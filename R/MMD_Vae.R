#' Compute MMD Penalty Using the IMQ Kernel
#'
#' This function computes the Maximum Mean Discrepancy (MMD) penalty between two tensors using the Inverse Multiquadric (IMQ) kernel.
#' It implements the unbiased U-statistic estimator of MMD with an IMQ kernel, as described in the Wasserstein Autoencoder (WAE)
#' framework. The IMQ kernel aggregates various scales by summing kernels computed at different scales, leveraging the property
#' that the sum of positive definite kernels remains positive definite.
#'
#' @param pz `tensorflow.tensor` First input tensor, typically representing samples from the model's latent space.
#' @param qz `tensorflow.tensor` Second input tensor, typically representing samples from a prior distribution.
#' @param batchSize `integer` Batch size used for computations (default is 32L).
#' @param sigmaZ `numeric` Kernel scaling factor (default is 1.0).
#' @param zDim `integer` Dimensionality of the latent space (default is 16L).
#' @return `tensorflow.tensor` A tensor containing the computed MMD penalty value.
#' @details This function is adapted from the Wasserstein Autoencoder (WAE) GitHub repository:
#'   <https://github.com/tolstikhin/wae/blob/master/wae.py#L233>. The penalty is calculated
#'   by summing contributions from IMQ kernels computed at multiple scales, enabling the model
#'   to capture discrepancies at various resolutions.
#' @examples
#' \dontrun{
#' library(tensorflow)
#' pz <- tf$random$normal(shape = c(32L, 16L))
#' qz <- tf$random$normal(shape = c(32L, 16L))
#' result <- mmdImqPenalty(pz, qz)
#' print(result)
#' }
#' @importFrom tensorflow tf
#' @export
mmdImqPenalty <- function(pz, qz, batchSize = 32L, sigmaZ = 1.0, zDim = 16L) {
  # Ensure input tensors have the correct shape
  if (!inherits(pz, "tensorflow.tensor") || !inherits(qz, "tensorflow.tensor")) {
    stop("Both 'pz' and 'qz' must be TensorFlow tensors.")
  }

  # Batch size inferred from tensor shape
  batchSize <- tf$shape(pz)[1]

  # Constants for computation
  nFloat <- tf$cast(batchSize, dtype = tf$float32)
  cBase <- 2 * zDim * sigmaZ

  # Compute pairwise squared distances for pz and qz
  computeDistances <- function(x, y) {
    normsX <- tf$reduce_sum(tf$square(x), axis = 1L, keepdims = TRUE)
    normsY <- tf$reduce_sum(tf$square(y), axis = 1L, keepdims = TRUE)
    dotProducts <- tf$matmul(x, y, transpose_b = TRUE)
    distances <- normsX + tf$transpose(normsY) - 2 * dotProducts
    return(distances)
  }

  distancesPz <- computeDistances(pz, pz)
  distancesQz <- computeDistances(qz, qz)
  distancesPq <- computeDistances(pz, qz)

  # Initialize MMD statistic
  mmdStat <- tf$constant(0.0, dtype = tf$float32)

  # Scales for the IMQ kernel
  scales <- c(0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0)

  # Compute MMD for each scale
  for (scale in scales) {
    C <- cBase * scale

    # Compute terms for the MMD statistic
    resPz <- C / (C + distancesPz)
    resPz <- resPz * (1 - tf$eye(batchSize, dtype = tf$float32))
    term1 <- tf$reduce_sum(resPz) / (nFloat * nFloat - nFloat)

    resQz <- C / (C + distancesQz)
    resQz <- resQz * (1 - tf$eye(batchSize, dtype = tf$float32))
    term2 <- tf$reduce_sum(resQz) / (nFloat * nFloat - nFloat)

    resPq <- C / (C + distancesPq)
    term3 <- 2 * tf$reduce_sum(resPq) / (nFloat * nFloat)

    # Update MMD statistic
    mmdStat <- mmdStat + (term1 + term2 - term3)
  }

  return(mmdStat)
}

#' Compute RBF Kernel for MMD Calculation
#'
#' Computes the Radial Basis Function (RBF) kernel between two input tensors, typically for use in Maximum Mean Discrepancy (MMD) calculations.
#'
#' @param x `tensorflow.tensor` First input tensor.
#' @param y `tensorflow.tensor` Second input tensor.
#' @return `tensorflow.tensor` A tensor containing the RBF kernel matrix.
#' @details This function calculates the RBF kernel by computing the element-wise squared differences between reshaped and tiled versions of the input tensors.
#' The result is scaled and exponentiated to form the kernel matrix.
#' @examples
#' \dontrun{
#' library(tensorflow)
#' x <- tf$random$normal(shape = c(32L, 16L))
#' y <- tf$random$normal(shape = c(32L, 16L))
#' kernelMatrix <- computeKernel(x, y)
#' print(kernelMatrix)
#' }
#' @importFrom tensorflow tf
#' @export
computeKernel <- function(x, y) {
  # Validate inputs
  if (!inherits(x, "tensorflow.tensor") || !inherits(y, "tensorflow.tensor")) {
    stop("Both 'x' and 'y' must be TensorFlow tensors.")
  }

  # Get tensor shapes
  xSize <- tf$shape(x)[1]
  ySize <- tf$shape(y)[1]
  dimensions <- tf$shape(x)[2]

  # Reshape and tile inputs for kernel computation
  tiledX <- tf$tile(tf$reshape(x, tf$stack(list(xSize, 1L, dimensions))), tf$stack(list(1L, ySize, 1L)))
  tiledY <- tf$tile(tf$reshape(y, tf$stack(list(1L, ySize, dimensions))), tf$stack(list(xSize, 1L, 1L)))

  # Compute RBF kernel
  kernelMatrix <- tf$exp(-tf$reduce_mean(tf$square(tiledX - tiledY), axis = 2L) / tf$cast(dimensions, tf$float32))

  return(kernelMatrix)
}

#' Compute Maximum Mean Discrepancy (MMD) Loss
#'
#' Calculates the MMD loss between two distributions using an RBF kernel.
#'
#' @param inputX `tensorflow.tensor` First input distribution.
#' @param inputY `tensorflow.tensor` Second input distribution.
#' @param sigmaSqr `numeric` Variance parameter for the RBF kernel, default is 1.0.
#' @return `tensorflow.tensor` Scalar MMD loss.
#' @examples
#' \dontrun{
#' library(tensorflow)
#' x <- tf$random$normal(shape = c(32L, 16L))
#' y <- tf$random$normal(shape = c(32L, 16L))
#' loss <- computeMMD(x, y)
#' print(loss)
#' }
#' @importFrom tensorflow tf
#' @export
computeMMD <- function(inputX, inputY, sigmaSqr = 1.0) {
  xKernel <- computeKernel(inputX, inputX)
  yKernel <- computeKernel(inputY, inputY)
  xyKernel <- computeKernel(inputX, inputY)

  tf$reduce_mean(xKernel) + tf$reduce_mean(yKernel) - 2 * tf$reduce_mean(xyKernel)
}

#' Train a Variational Autoencoder (VAE) with MMD Regularization
#'
#' Trains a VAE with MMD regularization, where the MMD loss enforces distributional
#' similarity between true and encoded samples in latent space.
#'
#' @param trainData `data.frame` Training data.
#' @param useMarkers `character` Names of marker columns to use.
#' @param epochs `integer` Number of training epochs, default is 80.
#' @param latentDim `integer` Dimensionality of latent space, optional.
#' @param seed `integer` Random seed for reproducibility, default is 1994.
#' @param lambdaWeight `numeric` Regularization weight for MMD loss, default is 0.1.
#' @param valData `data.frame` Validation data.
#' @param originalDim `integer` Dimensionality of input, default is 27.
#' @param batchSize `integer` Batch size for training, default is 16.
#' @param hiddenSizes `list` Optional hidden layer sizes.
#' @return `list` Containing the trained VAE model and encoder.
#' @examples
#' \dontrun{
#' trainedVAE <- trainVAEModel(trainData, useMarkers, valData = validationData)
#' }
#' @importFrom tensorflow tf set_random_seed
#' @importFrom keras layer_input layer_dense keras_model new_model_class
#' @export
trainVAEModel <- function(trainData, useMarkers, epochs = 80, latentDim = NULL, seed = 1994,
                          lambdaWeight = 0.1, valData, originalDim = 27L, batchSize = 16,
                          hiddenSizes = NULL) {
  tensorflow::set_random_seed(seed = seed)

  # Normalize and prepare data
  xTrain <- as.matrix(trainData[, useMarkers])
  xVal <- as.matrix(valData[, useMarkers])
  originalDim <- ncol(trainData)

  # Determine hidden layer sizes
  if (!is.null(hiddenSizes)) {
    intermediateDim <- hiddenSizes[[1]]
    intermediateDim2 <- hiddenSizes[[2]]
  } else {
    intermediateDim <- originalDim - 4
    intermediateDim2 <- intermediateDim - 4
  }

  if (is.null(latentDim)) latentDim <- intermediateDim2 - 3

  # Define encoder model
  encoderInputs <- layer_input(shape = originalDim)
  x <- encoderInputs %>%
    layer_dense(intermediateDim, activation = "relu") %>%
    layer_dense(intermediateDim2, activation = "relu") %>%
    layer_dense(latentDim, activation = "gelu")
  encoder <- keras_model(encoderInputs, x, name = "encoder")

  # Define decoder model
  decoderInputs <- layer_input(shape = latentDim)
  decoderOutputs <- decoderInputs %>%
    layer_dense(intermediateDim2, activation = "relu") %>%
    layer_dense(intermediateDim, activation = "relu") %>%
    layer_dense(originalDim, activation = "sigmoid")
  decoder <- keras_model(decoderInputs, decoderOutputs, name = "decoder")

  # Define custom VAE model
  modelVAE <- new_model_class(
    classname = "VAE",
    initialize = function(encoder, decoder, ...) {
      super$initialize(...)
      self$encoder <- encoder
      self$decoder <- decoder
      self$totalLossTracker <- metric_mean(name = "total_loss")
      self$reconstructionLossTracker <- metric_mean(name = "reconstruction_loss")
      self$mmdLossTracker <- metric_mean(name = "mmd_loss")
    },
    metrics = mark_active(function() {
      list(
        self$totalLossTracker,
        self$reconstructionLossTracker,
        self$mmdLossTracker
      )
    }),

    # Custom training step
    train_step = function(data) {
      x <- data[[1]]
      with(tf$GradientTape() %as% tape, {
        zMean <- self$encoder(x)
        reconstruction <- self$decoder(zMean)
        reconstructionLoss <- loss_binary_crossentropy(x, reconstruction) %>%
          op_sum(axis = 1) %>%
          op_mean()

        trueSamples <- tf$random$normal(shape = tf$shape(zMean))
        mmdLoss <- computeMMD(trueSamples, zMean)
        totalLoss <- reconstructionLoss + lambdaWeight * mmdLoss
      })

      grads <- tape$gradient(totalLoss, self$trainable_weights)
      self$optimizer$apply_gradients(zip_lists(grads, self$trainable_weights))

      self$totalLossTracker$update_state(totalLoss)
      self$reconstructionLossTracker$update_state(reconstructionLoss)
      self$mmdLossTracker$update_state(mmdLoss)

      list(
        total_loss = self$totalLossTracker$result(),
        reconstruction_loss = self$reconstructionLossTracker$result(),
        mmd_loss = self$mmdLossTracker$result()
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
      mmdLoss <- computeMMD(trueSamples, zMean)
      totalLoss <- reconstructionLoss + lambdaWeight * mmdLoss

      self$totalLossTracker$update_state(totalLoss)
      self$reconstructionLossTracker$update_state(reconstructionLoss)
      self$mmdLossTracker$update_state(mmdLoss)

      list(
        total_loss = self$totalLossTracker$result(),
        reconstruction_loss = self$reconstructionLossTracker$result(),
        mmd_loss = self$mmdLossTracker$result()
      )
    }
  )

  # Compile and train model
  vae <- modelVAE(encoder, decoder)

  lrSchedule <- learning_rate_schedule_exponential_decay(
    1e-3,
    decay_steps = 100000, decay_rate = 0.95, staircase = FALSE
  )
  opt <- optimizer_rmsprop(learning_rate = lrSchedule, momentum = 0, centered = TRUE)
  vae %>% compile(optimizer = opt)

  esCallback <- callback_early_stopping(
    min_delta = 1e-4, monitor = "val_total_loss", mode = "min",
    patience = 15, verbose = 1, restore_best_weights = TRUE
  )

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
#' Encodes and decodes new samples through the Variational Autoencoder (VAE) model, returning both
#' latent representations and decoded samples.
#'
#' @param newSamples `data.frame` New input samples to encode and decode.
#' @param vae `list` Trained VAE model object containing the encoder and decoder.
#' @param latentDim `integer` Dimensionality of the latent space, default is 8.
#' @param batchSize `integer` Batch size for decoding, default is 16.
#' @return `list` A list with two components:
#'   - `decoded`: Decoded samples as a data frame.
#'   - `encoded`: Latent representations as a tensor.
#' @examples
#' \dontrun{
#' # Assuming `vaeModel` is a trained VAE and `testSamples` is a data frame of new samples:
#' result <- decodeSamples(newSamples = testSamples, vae = vaeModel)
#' print(result$decoded) # Decoded samples
#' print(result$encoded) # Latent representations
#' }
#' @importFrom tensorflow tf
#' @importFrom keras predict
#' @export
decodeSamples <- function(newSamples, vae, latentDim = 8L, batchSize = 16) {
  # Set random seed for reproducibility
  tensorflow::set_random_seed(seed = 1994)

  # Validate inputs
  if (!inherits(newSamples, "data.frame")) {
    stop("'newSamples' must be a data frame.")
  }
  if (!is.list(vae) || is.null(vae$encoder) || is.null(vae$decoder)) {
    stop("'vae' must be a list containing 'encoder' and 'decoder' models.")
  }

  # Encode the input samples
  zMean <- predict(vae$encoder, newSamples, batch_size = batchSize)

  # Decode the latent representations
  decodedSamples <- predict(vae$decoder, zMean, batch_size = batchSize) %>%
    as.data.frame()

  # Return encoded and decoded results
  list(decoded = decodedSamples, encoded = zMean)
}
