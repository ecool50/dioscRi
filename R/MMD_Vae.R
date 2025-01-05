    mmd_penalty <- function(pz, qz, batch_size = 32L,
                            sigma_z = 1., z_dim = 16L) {
      # This method calculates the unbiased U-statistic estimator of
      #         the MMD with the IMQ kernel. It's taken from
      #         https://github.com/tolstikhin/wae/blob/master/wae.py#L233
      #
      #         Here the property that the sum of positive definite kernels is
      #         still a p.d. kernel is used. Various kernels calculated at different
      #         scales are summed together in order to "simultaneously look at various
      #   scales [https://github.com/tolstikhin/wae/issues/2].

      batch_size <- tf$shape(pz)[1]
      nf <- tf$cast(batch_size, dtype = tf$float32)


      norms_pz <- tf$reduce_sum(tf$square(pz), axis = 1L, keepdims = TRUE)
      dotprods_pz <- tf$matmul(pz, qz, transpose_b = TRUE)
      distances_pz <- norms_pz + tf$transpose(norms_pz) - 2. * dotprods_pz

      norms_qz <- tf$reduce_sum(tf$square(qz), axis = 1L, keepdims = TRUE)
      dotprods_qz <- tf$matmul(qz, qz, transpose_b = TRUE)
      distances_qz <- norms_qz + tf$transpose(norms_qz) - 2. * dotprods_qz

      dotprods <- tf$matmul(qz, pz, transpose_b = TRUE)
      distances <- norms_qz + tf$transpose(norms_pz) - 2. * dotprods

      cbase <- tf$constant(2. * z_dim * sigma_z)
      stat <- tf$constant(0.)
      nf <- tf$cast(batch_size, dtype = tf$float32)

      scales <- c(0.1, 0.2, 0.5, 1., 2., 5., 10.)
      for (scale in scales) {
        C <- cbase * scale
        res1 <- C / (C + distances_qz)
        res1 <- res1 + (C / (C + distances_pz))
        res1 <- tf$multiply(res1, 1. - tf$eye(batch_size, dtype = tf$float32))
        res1 <- tf$reduce_sum(res1) / (nf * nf - nf)

        res2 <- C / (C + distances)
        res2 <- tf$reduce_sum(res2) * 2. / (nf * nf)

        stat <- stat + (res1 - res2)
      }

      stat
    }

    #' Compute RBF Kernel for MMD Calculation
    #'
    #' Computes the RBF kernel between two inputs, typically for use in MMD calculations.
    #' @param x Tensor; first input tensor.
    #' @param y Tensor; second input tensor.
    #' @return Tensor containing the RBF kernel matrix.
    #' @importFrom tensorflow tf
    #' @noRd
    computeKernel <- function(x, y) {
      xSize <- tf$shape(x)[1]
      ySize <- tf$shape(y)[1]
      dim <- tf$shape(x)[2]

      # Reshape and tile inputs for kernel computation
      tiledX <- tf$tile(tf$reshape(x, tf$stack(list(xSize, 1L, dim))), tf$stack(list(1L, ySize, 1L)))
      tiledY <- tf$tile(tf$reshape(y, tf$stack(list(1L, ySize, dim))), tf$stack(list(xSize, 1L, 1L)))

      tf$exp(-tf$reduce_mean(tf$square(tiledX - tiledY), axis = 2L) / tf$cast(dim, tf$float32))
    }

    #' Compute Maximum Mean Discrepancy (MMD) Loss
    #'
    #' Calculates the MMD loss between two distributions using an RBF kernel.
    #' @param x Tensor; first input distribution.
    #' @param y Tensor; second input distribution.
    #' @param sigmaSqr Numeric; variance parameter for RBF kernel, default is 1.0.
    #' @return Tensor; scalar MMD loss.
    #' @importFrom tensorflow tf
    #' @noRd
    computeMMD <- function(x, y, sigmaSqr = 1.0) {
      xKernel <- computeKernel(x, x)
      yKernel <- computeKernel(y, y)
      xyKernel <- computeKernel(x, y)

      tf$reduce_mean(xKernel) + tf$reduce_mean(yKernel) - 2 * tf$reduce_mean(xyKernel)
    }

    #' Train a Variational Autoencoder (VAE) with MMD Regularization
    #'
    #' Trains a VAE with MMD regularization, where the MMD loss enforces distributional
    #' similarity between true and encoded samples in latent space.
    #' @param trainData Data frame of training data.
    #' @param useMarkers Character vector; names of marker columns to use.
    #' @param epochs Integer; number of training epochs, default is 80.
    #' @param latentDim Integer; dimensionality of latent space, optional.
    #' @param lambda Numeric; regularization weight for MMD loss, default is 0.1.
    #' @param valData Data frame of validation data.
    #' @param originalDim Integer; dimensionality of input, default is 27.
    #' @param batchSize Integer; batch size for training, default is 16.
    #' @return List containing the trained VAE model and encoder.
    #' @importFrom tensorflow tf set_random_seed
    #' @importFrom keras layer_input layer_dense keras_model new_model_class
    #' @noRd
    trainVAEModel <- function(trainData, useMarkers, epochs = 80, latentDim = NULL, seed = 1994,
                              lambda = 0.1, valData, originalDim = 27L, batchSize = 16,
                              hiddenSizes = NULL) {
      tensorflow::set_random_seed(seed = seed)

      # Normalize and prepare data
      xTrain <- as.matrix(trainData[, useMarkers])
      xVal <- as.matrix(valData[, useMarkers])
      originalDim <- ncol(trainData)

      if (!is.null(hiddenSizes)) {
        intermediateDim <- hiddenSizes[[1]]
        intermediateDim2 <- hiddenSizes[[2]]
      } else {
        # Define model dimensions
        intermediateDim <- originalDim - 4
        intermediateDim2 <- intermediateDim - 4
      }


      if (is.null(latentDim)) latentDim <- intermediateDim2 - 3

      # Encoder model
      encoderInputs <- layer_input(shape = originalDim)
      x <- encoderInputs %>%
        layer_dense(intermediateDim, activation = "relu") %>%
        layer_dense(intermediateDim2, activation = "relu") %>%
        layer_dense(latentDim, activation = "gelu")
      encoder <- keras_model(encoderInputs, x, name = "encoder")

      # Decoder model
      decoderInputs <- layer_input(shape = latentDim)
      decoderOutputs <- decoderInputs %>%
        layer_dense(intermediateDim2, activation = "relu") %>%
        layer_dense(intermediateDim, activation = "relu") %>%
        layer_dense(originalDim, activation = "sigmoid")
      decoder <- keras_model(decoderInputs, decoderOutputs, name = "decoder")

      # Custom VAE model class
      modelVAE <- new_model_class(
        classname = "VAE",
        initialize = function(encoder, decoder, ...) {
          super$initialize(...)
          self$encoder <- encoder
          self$decoder <- decoder
          self$total_loss_tracker <- metric_mean(name = "total_loss")
          self$reconstruction_loss_tracker <- metric_mean(name = "reconstruction_loss")
          self$mmd_loss_tracker <- metric_mean(name = "mmd_loss")
        },
        metrics = mark_active(function() {
          list(
            self$total_loss_tracker,
            self$reconstruction_loss_tracker,
            self$mmd_loss_tracker
          )
        }),

        # Custom training step with MMD regularization
        train_step = function(data) {
          x <- data[[1]]
          with(tf$GradientTape() %as% tape, {
            zMean <- self$encoder(x)
            reconstruction <- self$decoder(zMean)
            reconstructionLoss <- loss_binary_crossentropy(x, reconstruction) %>%
              op_sum(axis = 1) %>%
              op_mean()

            trueSamples <- tf$random$normal(shape = tf$shape(zMean))
            mmdLoss <- mmd_penalty(trueSamples, zMean)
            totalLoss <- reconstructionLoss + lambda * mmdLoss
          })

          grads <- tape$gradient(totalLoss, self$trainable_weights)
          self$optimizer$apply_gradients(zip_lists(grads, self$trainable_weights))

          self$total_loss_tracker$update_state(totalLoss)
          self$reconstruction_loss_tracker$update_state(reconstructionLoss)
          self$mmd_loss_tracker$update_state(mmdLoss)

          list(
            total_loss = self$total_loss_tracker$result(),
            reconstruction_loss = self$reconstruction_loss_tracker$result(),
            mmd_loss = self$mmd_loss_tracker$result()
          )
        },

        # Custom validation step
        test_step = function(data) {
          x <- data[[1]]
          zMean <- self$encoder(x)
          reconstruction <- self$decoder(zMean)
          reconstructionLoss <- loss_binary_crossentropy(x, reconstruction) %>%
            op_sum(axis = 1) %>%
            op_mean()

          trueSamples <- tf$random$normal(shape = tf$shape(zMean))
          mmdLoss <- mmd_penalty(trueSamples, zMean)
          totalLoss <- reconstructionLoss + lambda * mmdLoss

          self$total_loss_tracker$update_state(totalLoss)
          self$reconstruction_loss_tracker$update_state(reconstructionLoss)
          self$mmd_loss_tracker$update_state(mmdLoss)

          list(
            total_loss = self$total_loss_tracker$result(),
            reconstruction_loss = self$reconstruction_loss_tracker$result(),
            mmd_loss = self$mmd_loss_tracker$result()
          )
        }
      )

      # Instantiate and compile VAE model
      vae <- modelVAE(encoder, decoder)

      # Define learning rate schedule and optimizer
      lrSchedule <- learning_rate_schedule_exponential_decay(
        1e-3,
        decay_steps = 100000, decay_rate = 0.95, staircase = FALSE
      )
      opt <- optimizer_rmsprop(learning_rate = lrSchedule, momentum = 0, centered = TRUE)
      vae %>% compile(optimizer = opt)

      # Early stopping callback
      esCallback <- callback_early_stopping(
        min_delta = 1e-4, monitor = "val_total_loss", mode = "min",
        patience = 15, verbose = 1, restore_best_weights = TRUE
      )

      # Fit model with validation data
      vae %>% fit(
        xTrain, xTrain,
        batch_size = batchSize,
        epochs = epochs,
        validation_data = list(xVal, xVal),
        shuffle = TRUE
      )

      return(list(vae = vae, encoder = encoder))
    }

    #' Decode New Samples Using VAE Decoder
    #'
    #' Encodes and decodes new samples through the VAE model, returning both
    #' latent representations and decoded samples.
    #' @param newSamples Data frame of new input samples to encode and decode.
    #' @param vae Trained VAE model object.
    #' @param latentDim Integer; dimensionality of the latent space, default is 8.
    #' @param batchSize Integer; batch size for decoding, default is 16.
    #' @return List with decoded samples and encoded latent representations.
    #' @importFrom tensorflow tf
    #' @importFrom keras predict
    #' @noRd
    decodeSamples <- function(newSamples, vae, latentDim = 8L, batchSize = 16) {
      tensorflow::set_random_seed(seed = 1994)
      zMean <- predict(vae$encoder, newSamples)

      # Decode the latent representation
      decodedSamples <- predict(vae$decoder, zMean) %>% as.data.frame()

      list(decoded = decodedSamples, encoded = zMean)
    }
