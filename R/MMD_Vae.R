#' Compute Maximum Mean Discrepancy (MMD) Penalty with IMQ Kernel
#'
#' This function calculates the unbiased U-statistic estimator of the Maximum Mean Discrepancy (MMD)
#' using the Inverse Multiquadric (IMQ) kernel. The IMQ kernel aggregates various scales by summing 
#' kernels computed at different scales, leveraging the property that the sum of positive definite 
#' kernels remains positive definite. This approach enables the model to analyze discrepancies 
#' across multiple resolutions.
#'
#' @param pz `tensorflow.tensor` First input tensor, typically representing encoded samples.
#' @param qz `tensorflow.tensor` Second input tensor, typically representing samples from a prior distribution.
#' @param batchSize `integer` Batch size used for computation, default is 32L.
#' @param sigmaZ `numeric` Scaling factor for the IMQ kernel, default is 1.0.
#' @param zDim `integer` Dimensionality of the latent space, default is 16L.
#' @return `tensorflow.tensor` A scalar tensor containing the computed MMD penalty value.
#' @details This implementation is adapted from the Wasserstein Autoencoder (WAE) repository:
#'   <https://github.com/tolstikhin/wae/blob/master/wae.py#L233>. The penalty is calculated 
#'   by summing contributions from IMQ kernels computed at multiple scales.
#' @importFrom tensorflow tf
#' @export
#' @noRd
mmdPenalty <- function(pz, qz, batch_size = 32L,
                        sigmaZ = 1., zDim = 16L){
    
    # This method calculates the unbiased U-statistic estimator of
    #         the MMD with the IMQ kernel. It's taken from
    #         https://github.com/tolstikhin/wae/blob/master/wae.py#L233
    #
    #         Here the property that the sum of positive definite kernels is
    #         still a p.d. kernel is used. Various kernels calculated at different
    #         scales are summed together in order to "simultaneously look at various
    #   scales [https://github.com/tolstikhin/wae/issues/2].
    
    batchSize <- tf$shape(pz)[1]
    nf <- tf$cast(batchSize, dtype = tf$float32)
    
    
    normsPz <- tf$reduce_sum(tf$square(pz), axis = 1L, keepdims=TRUE)
    dotProdsPz <- tf$matmul(pz,qz, transpose_b=TRUE)
    distancesPz <- normsPz + tf$transpose(normsPz) - 2. * dotProdsPz
    
    normsQz <- tf$reduce_sum(tf$square(qz), axis = 1L, keepdims=TRUE)
    dotProdsQz <- tf$matmul(qz,qz, transpose_b=TRUE)
    distancesQz <- normsQz + tf$transpose(normsQz) - 2. * dotProdsQz
    
    dotProds <- tf$matmul(qz, pz, transpose_b=TRUE)
    distances <- normsQz + tf$transpose(normsPz) - 2. * dotProds
    
    cBase <- tf$constant(2. * zDim * sigmaZ)
    stat <- tf$constant(0.)
    nf <- tf$cast(batchSize, dtype=tf$float32)
    
    scales <- c(0.1, 0.2, 0.5, 1., 2., 5., 10.)
    for (scale in scales) {
        C <- cBase * scale
        res1 <- C / (C + distancesQz)
        res1 <- res1 + (C / (C + distancesPz))
        res1 <- tf$multiply(res1, 1. - tf$eye(batchSize, dtype = tf$float32))
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
#' @export
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
#' @export
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
#' @importFrom keras3 layer_input layer_dense keras_model new_model_class
#' @export
#' @examples
#' # Load sample data
#' data(sample_cytof_data)
#' data(sample_markers)
#' 
#' # First compute reference samples for optimal training/validation split
#' ref_samples <- computeReferenceSample(
#'   data = sample_cytof_data,
#'   markers = sample_markers,
#'   sampleCol = "sample_id",
#'   N = 2
#' )
#' 
#' # Use top samples (most representative) for training
#' train_data <- sample_cytof_data[
#'   sample_cytof_data$sample_id %in% ref_samples$topNSamples,
#'   sample_markers
#' ]
#' 
#' # Use bottom samples (most distinct) for validation
#' val_data <- sample_cytof_data[
#'   sample_cytof_data$sample_id %in% ref_samples$bottomNSamples,
#'   sample_markers
#' ]
#' 
#' # Train VAE model with optimal data split
#' vae_model <- trainVAEModel(
#'   trainData = train_data,
#'   useMarkers = sample_markers,
#'   epochs = 30,
#'   latentDim = 2,
#'   hiddenSizes = c(4, 3),
#'   lambda = 0.01,
#'   valData = val_data,
#'   originalDim = length(sample_markers),
#'   batchSize = 32,
#'   seed = 1994
#' )
#' 
#' # Access components
#' encoder <- vae_model$encoder
#' vae <- vae_model$vae
trainVAEModel <- function(trainData, useMarkers, epochs = 80, latentDim = NULL, seed = 1994,
                          lambda = 0.1, valData, originalDim = 27L, batchSize = 16,
                          hiddenSizes = NULL) {
    
    tensorflow::set_random_seed(seed = seed)
    
    # Normalize and prepare data
    xTrain <- as.matrix(trainData[, useMarkers])
    xVal <- as.matrix(valData[, useMarkers])
    originalDim <- ncol(trainData)
    
    if(!is.null(hiddenSizes)){
        intermediateDim <- hiddenSizes[[1]]
        intermediateDim2 <- hiddenSizes[[2]]
    }else{
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
                mmdLoss <- mmdPenalty(trueSamples, zMean)
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
            mmdLoss <- mmdPenalty(trueSamples, zMean)
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
        1e-3, decay_steps = 100000, decay_rate = 0.95, staircase = FALSE
    )
    opt <- optimizer_rmsprop(learning_rate = lrSchedule, momentum = 0, centered = TRUE)
    vae %>% compile(optimizer = opt)
    
    # Early stopping callback
    esCallback <- callback_early_stopping(
        min_delta = 1e-4, monitor = 'val_total_loss', mode = 'min',
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
#' @export
#' @examples
#' # Assuming vae_model is already trained (from previous example)
#' data(sample_cytof_data)
#' data(sample_markers)
#' 
#' # Normalize new samples
#' new_data <- sample_cytof_data[1:100, sample_markers]
#' 
#' normalized <- decodeSamples(
#'   newSamples = as.matrix(new_data),
#'   vae = vae_model$vae,
#'   latentDim = 5,
#'   batchSize = 32
#' )
#' 
#' # Access results
#' decoded_data <- normalized$decoded
#' latent_rep <- normalized$encoded
#' 
#' dim(decoded_data)
#' dim(latent_rep)
decodeSamples <- function(newSamples, vae, latentDim = 16L, batchSize = 32L) {
    tensorflow::set_random_seed(seed = 1994)
    zMean <- predict(vae$encoder, newSamples)
    
    # Decode the latent representation
    decodedSamples <- predict(vae$decoder, zMean) %>% as.data.frame()
    
    list(decoded = decodedSamples, encoded = zMean)
}


#' Compute Sliced Wasserstein Distance
#'
#' Computes the Sliced Wasserstein Distance between two distributions by projecting
#' onto random 1D directions and computing the 1D Wasserstein distance (sorted L2).
#'
#' @param z `tensorflow.tensor` Encoded samples from the autoencoder.
#' @param priorZ `tensorflow.tensor` Samples from the prior distribution (e.g., N(0,1)).
#' @param numProjections `integer` Number of random projection directions, default is 50L.
#' @return `tensorflow.tensor` A scalar tensor containing the sliced Wasserstein distance.
#' @importFrom tensorflow tf
#' @noRd
slicedWassersteinDistance <- function(z, priorZ, numProjections = 50L) {
    tf <- tensorflow::tf

    # Get dimensions
    latentDim <- tf$shape(z)[2]
    batchSize <- tf$shape(z)[1]

    # Ensure numProjections is integer
    numProjections <- as.integer(numProjections)

    # Generate random projection directions (unit vectors on hypersphere)
    # Shape: [latentDim, numProjections]
    theta <- tf$random$normal(shape = list(latentDim, numProjections))
    theta <- theta / tf$norm(theta, axis = 0L, keepdims = TRUE)

    # Project both distributions onto random directions
    # z: [batch, latentDim] @ theta: [latentDim, numProjections] -> [batch, numProjections]
    zProj <- tf$matmul(z, theta)
    priorProj <- tf$matmul(priorZ, theta)

    # Sort along batch dimension for each projection
    zSorted <- tf$sort(zProj, axis = 0L)
    priorSorted <- tf$sort(priorProj, axis = 0L)

    # Compute L2 Sliced Wasserstein Distance
    # Mean over all projections and all sorted pairs
    tf$reduce_mean(tf$square(zSorted - priorSorted))
}


#' Compute Distributional Sliced Wasserstein Distance
#'
#' Computes the Distributional Sliced Wasserstein Distance (DSW) between two
#' distributions. Unlike vanilla SW which uses uniform random projections, DSW
#' learns an optimal distribution over projections that focuses on the most
#' discriminative directions. Based on Nguyen et al. (2021).
#'
#' @param z `tensorflow.tensor` Encoded samples from the autoencoder.
#' @param priorZ `tensorflow.tensor` Samples from the prior distribution.
#' @param numProjections `integer` Number of projection directions.
#' @param projectionNet Keras model that transforms random projections to optimized ones.
#' @param projOptimizer Optimizer for the projection network.
#' @param maxIter `integer` Number of inner optimization steps, default 10.
#' @param lam `numeric` Regularization strength for projection diversity, default 1.0.
#' @return `tensorflow.tensor` A scalar tensor containing the DSW distance.
#' @importFrom tensorflow tf
#' @noRd
distributionalSlicedWassersteinDistance <- function(z, priorZ, numProjections = 50L,
                                                    projectionNet, projOptimizer,
                                                    maxIter = 10L, lam = 1.0) {
    tf <- tensorflow::tf

    latentDim <- tf$shape(z)[2]
    numProjections <- as.integer(numProjections)
    maxIter <- as.integer(maxIter)

    # Generate initial random projections on unit sphere
    pro <- tf$random$normal(shape = list(numProjections, latentDim))
    pro <- pro / tf$norm(pro, axis = 1L, keepdims = TRUE)

    # Detach samples for inner optimization (don't backprop to encoder)
    zDetach <- tf$stop_gradient(z)
    priorDetach <- tf$stop_gradient(priorZ)

    # Inner optimization: find projections that maximize SW distance
    for (i in seq_len(maxIter)) {
        with(tf$GradientTape() %as% tape, {
            # Transform random projections through learned network
            projections <- projectionNet(pro, training = TRUE)
            projections <- projections / tf$norm(projections, axis = 1L, keepdims = TRUE)

            # Cosine similarity regularization: encourage diverse projections
            cosine <- tf$matmul(projections, projections, transpose_b = TRUE)
            identity <- tf$eye(numProjections)
            reg <- lam * tf$reduce_mean(tf$square(cosine - identity))

            # Compute 1D Wasserstein distances
            # z: [batch, latent], projections: [numProj, latent]
            zProj <- tf$matmul(zDetach, projections, transpose_b = TRUE)
            priorProj <- tf$matmul(priorDetach, projections, transpose_b = TRUE)

            zSorted <- tf$sort(zProj, axis = 0L)
            priorSorted <- tf$sort(priorProj, axis = 0L)

            wd <- tf$reduce_mean(tf$square(zSorted - priorSorted))

            # Maximize WD while keeping projections diverse
            innerLoss <- reg - wd
        })

        grads <- tape$gradient(innerLoss, projectionNet$trainable_variables)
        projOptimizer$apply_gradients(Map(list, grads, projectionNet$trainable_variables))
    }

    # Final computation with gradients flowing to encoder
    projections <- projectionNet(pro, training = FALSE)
    projections <- projections / tf$norm(projections, axis = 1L, keepdims = TRUE)

    zProj <- tf$matmul(z, projections, transpose_b = TRUE)
    priorProj <- tf$matmul(priorZ, projections, transpose_b = TRUE)

    zSorted <- tf$sort(zProj, axis = 0L)
    priorSorted <- tf$sort(priorProj, axis = 0L)

    tf$reduce_mean(tf$square(zSorted - priorSorted))
}


#' Train a Sliced Wasserstein Autoencoder (SWAE)
#'
#' Trains an autoencoder with Sliced Wasserstein Distance regularization.
#' SWAE uses optimal transport geometry to match the latent distribution to the prior,
#' offering computational advantages over MMD-based approaches (O(n log n) vs O(n^2)).
#'
#' @param trainData Data frame of training data.
#' @param useMarkers Character vector; names of marker columns to use.
#' @param epochs Integer; number of training epochs, default is 80.
#' @param latentDim Integer; dimensionality of latent space, optional.
#' @param lambda Numeric; regularization weight for Sliced Wasserstein loss, default is 10.0.
#' @param numProjections Integer; number of random projections for SW distance, default is 50L.
#' @param valData Data frame of validation data.
#' @param originalDim Integer; dimensionality of input, default is 27.
#' @param batchSize Integer; batch size for training, default is 32.
#' @param hiddenSizes Numeric vector; sizes of hidden layers in encoder/decoder.
#' @param seed Integer; random seed for reproducibility, default is 1994.
#' @return List containing the trained SWAE model, encoder, and decoder.
#' @importFrom tensorflow tf set_random_seed
#' @importFrom keras3 layer_input layer_dense keras_model
#' @export
#' @examples
#' \dontrun{
#' # Requires TensorFlow and keras3 to be installed
#' data(sample_cytof_data)
#' data(sample_markers)
#'
#' ref_samples <- computeReferenceSample(
#'   data = sample_cytof_data,
#'   markers = sample_markers,
#'   sampleCol = "sample_id",
#'   N = 2
#' )
#'
#' train_data <- sample_cytof_data[
#'   sample_cytof_data$sample_id %in% ref_samples$topNSamples,
#'   sample_markers
#' ]
#' val_data <- sample_cytof_data[
#'   sample_cytof_data$sample_id %in% ref_samples$bottomNSamples,
#'   sample_markers
#' ]
#'
#' swae_model <- trainSWAEModel(
#'   trainData = train_data,
#'   useMarkers = sample_markers,
#'   epochs = 30,
#'   latentDim = 8,
#'   lambda = 10.0,
#'   numProjections = 50,
#'   valData = val_data,
#'   batchSize = 32,
#'   seed = 1994
#' )
#'
#' # Use with decodeSamples (same as WAE-MMD)
#' normalized <- decodeSamples(newSamples = as.matrix(new_data), vae = swae_model)
#' }
trainSWAEModel <- function(trainData, useMarkers, epochs = 80, latentDim = NULL, seed = 1994,
                           lambda = 10.0, numProjections = 50L, valData, originalDim = 27L,
                           batchSize = 32, hiddenSizes = NULL) {

    # Set seed for reproducibility
    tensorflow::set_random_seed(seed = seed)
    tf <- tensorflow::tf

    # Prepare data
    xTrain <- as.matrix(trainData[, useMarkers])
    xVal <- as.matrix(valData[, useMarkers])
    originalDim <- ncol(xTrain)  # Use marker count, not all columns

    # Determine hidden layer sizes
    if (!is.null(hiddenSizes)) {
        intermediateDim <- hiddenSizes[[1]]
        intermediateDim2 <- hiddenSizes[[2]]
    } else {
        intermediateDim <- max(originalDim - 4, ceiling(originalDim * 0.75), 4)
        intermediateDim2 <- max(intermediateDim - 4, ceiling(intermediateDim * 0.75), 3)
    }

    if (is.null(latentDim)) latentDim <- max(intermediateDim2 - 3, 2)

    # Ensure integer types
    numProjections <- as.integer(numProjections)
    latentDim <- as.integer(latentDim)

    message(sprintf("SWAE Architecture: %d -> %d -> %d -> %d -> %d -> %d -> %d",
                    originalDim, intermediateDim, intermediateDim2, latentDim,
                    intermediateDim2, intermediateDim, originalDim))
    message(sprintf("Using %d random projections for Sliced Wasserstein Distance", numProjections))

    # Encoder model
    encoderInputs <- keras3::layer_input(shape = originalDim)
    x <- encoderInputs |>
        keras3::layer_dense(intermediateDim, activation = "relu") |>
        keras3::layer_dense(intermediateDim2, activation = "relu") |>
        keras3::layer_dense(latentDim, activation = "gelu")
    encoder <- keras3::keras_model(encoderInputs, x, name = "encoder")

    # Decoder model
    decoderInputs <- keras3::layer_input(shape = latentDim)
    decoderOutputs <- decoderInputs |>
        keras3::layer_dense(intermediateDim2, activation = "relu") |>
        keras3::layer_dense(intermediateDim, activation = "relu") |>
        keras3::layer_dense(originalDim, activation = "sigmoid")
    decoder <- keras3::keras_model(decoderInputs, decoderOutputs, name = "decoder")

    # Build functional model
    inputs <- keras3::layer_input(shape = originalDim)
    z <- encoder(inputs)
    outputs <- decoder(z)
    swae <- keras3::keras_model(inputs, outputs, name = "swae")

    # Custom loss function combining reconstruction and Sliced Wasserstein
    swae_loss <- function(y_true, y_pred) {
        # Reconstruction loss (MSE works well for SWAE)
        reconstruction_loss <- tf$reduce_mean(
            tf$square(y_true - y_pred)
        ) * tf$cast(originalDim, tf$float32)

        # Get latent representation
        z_encoded <- encoder(y_true)

        # Sample from prior N(0, 1)
        prior_samples <- tf$random$normal(shape = tf$shape(z_encoded))

        # Sliced Wasserstein distance
        sw_loss <- slicedWassersteinDistance(z_encoded, prior_samples, numProjections)

        # Total loss
        reconstruction_loss + lambda * sw_loss
    }

    # Compile model with Adam optimizer
    swae |> keras3::compile(
        optimizer = keras3::optimizer_adam(learning_rate = 1e-3),
        loss = swae_loss
    )

    # Early stopping callback
    esCallback <- keras3::callback_early_stopping(
        min_delta = 1e-4,
        monitor = "val_loss",
        mode = "min",
        patience = 15,
        verbose = 1,
        restore_best_weights = TRUE
    )

    # Learning rate reduction callback
    lrCallback <- keras3::callback_reduce_lr_on_plateau(
        monitor = "val_loss",
        factor = 0.5,
        patience = 5,
        min_lr = 1e-6,
        verbose = 1
    )

    # Fit model
    history <- swae |> keras3::fit(
        xTrain, xTrain,
        batch_size = batchSize,
        epochs = epochs,
        validation_data = list(xVal, xVal),
        shuffle = TRUE,
        callbacks = list(esCallback, lrCallback),
        verbose = 1
    )

    return(list(vae = swae, encoder = encoder, decoder = decoder, history = history))
}


#' Train Enhanced Sliced Wasserstein Autoencoder (SWAE) with Annealing
#'
#' Trains an autoencoder with Sliced Wasserstein Distance regularization using
#' a custom training loop that supports lambda annealing. This version offers
#' better control over training dynamics and typically produces better latent
#' representations than the standard SWAE.
#'
#' @param trainData Data frame or matrix of training data.
#' @param useMarkers Character vector; names of marker columns to use.
#' @param epochs Integer; number of training epochs, default is 100.
#' @param latentDim Integer; dimensionality of latent space. Default is 32.
#' @param lambda Numeric; final regularization weight for SW loss, default is 10.0.
#' @param lambdaAnnealEpochs Integer; epochs to linearly anneal lambda from 0 to final value.
#'   Default is 0 (no annealing). Recommended: 20-50 epochs for better representations.
#' @param numProjections Integer; number of random projections for SW distance, default is 100L.
#' @param batchSize Integer; batch size for training, default is 256.
#' @param learningRate Numeric; learning rate, default is 1e-3.
#' @param hiddenDims Integer vector; sizes of hidden layers. Default is c(256, 128).
#' @param valData Optional validation data for monitoring.
#' @param verbose Integer; print frequency. Default is 10 (every 10 epochs).
#' @param seed Integer; random seed for reproducibility, default is 1994.
#' @return List containing encoder, decoder, and training history.
#' @importFrom tensorflow tf set_random_seed
#' @importFrom keras3 layer_input layer_dense keras_model
#' @importFrom reticulate %as%
#' @export
#' @examples
#' \dontrun{
#' # Train enhanced SWAE with lambda annealing
#' swae_model <- trainSWAEEnhanced(
#'   trainData = train_data,
#'   useMarkers = markers,
#'   epochs = 100,
#'   latentDim = 32,
#'   lambda = 10.0,
#'   lambdaAnnealEpochs = 30,  # Anneal over first 30 epochs
#'   hiddenDims = c(256, 128),
#'   batchSize = 256,
#'   seed = 1994
#' )
#'
#' # Get latent representations
#' latent <- encodeSWAE(swae_model, new_data)
#' }
trainSWAEEnhanced <- function(trainData, useMarkers, epochs = 100L, latentDim = 32L,
                               lambda = 10.0, lambdaAnnealEpochs = 0L, numProjections = 100L,
                               batchSize = 256L, learningRate = 1e-3,
                               hiddenDims = c(256L, 128L), valData = NULL,
                               verbose = 10L, seed = 1994) {

    tensorflow::set_random_seed(seed)
    tf <- tensorflow::tf
    set.seed(seed)

    # Prepare data
    xTrain <- as.matrix(trainData[, useMarkers])
    inputDim <- ncol(xTrain)
    nTrain <- nrow(xTrain)

    # Ensure integers
    latentDim <- as.integer(latentDim)
    batchSize <- as.integer(batchSize)
    numProjections <- as.integer(numProjections)
    hiddenDims <- as.integer(hiddenDims)
    lambdaAnnealEpochs <- as.integer(lambdaAnnealEpochs)

    message("Enhanced SWAE Architecture:")
    message(sprintf("  Input dim: %d", inputDim))
    message(sprintf("  Hidden dims: %s", paste(hiddenDims, collapse = " -> ")))
    message(sprintf("  Latent dim: %d", latentDim))
    message(sprintf("  Training cells: %d", nTrain))
    message(sprintf("  Lambda (SW): %.2f%s", lambda,
                   ifelse(lambdaAnnealEpochs > 0, sprintf(" (annealed over %d epochs)", lambdaAnnealEpochs), "")))
    message(sprintf("  Projections: %d", numProjections))

    # ========== BUILD ENCODER ==========
    encoderInput <- keras3::layer_input(shape = inputDim, name = "encoder_input")

    h <- encoderInput
    for (i in seq_along(hiddenDims)) {
        h <- h |> keras3::layer_dense(units = hiddenDims[i], activation = "relu",
                                       name = paste0("encoder_h", i))
    }
    # Linear output for latent space (SWAE doesn't need reparameterization)
    latentOutput <- h |> keras3::layer_dense(units = latentDim, activation = "gelu",
                                              name = "latent")

    encoder <- keras3::keras_model(encoderInput, latentOutput, name = "encoder")

    # ========== BUILD DECODER ==========
    decoderInput <- keras3::layer_input(shape = latentDim, name = "decoder_input")

    h <- decoderInput
    for (i in rev(seq_along(hiddenDims))) {
        h <- h |> keras3::layer_dense(units = hiddenDims[i], activation = "relu",
                                       name = paste0("decoder_h", length(hiddenDims) - i + 1))
    }
    decoderOutput <- h |> keras3::layer_dense(units = inputDim, activation = "sigmoid",
                                               name = "decoder_output")

    decoder <- keras3::keras_model(decoderInput, decoderOutput, name = "decoder")

    # ========== COLLECT TRAINABLE VARIABLES ==========
    allVariables <- c(encoder$trainable_variables, decoder$trainable_variables)

    # ========== OPTIMIZER ==========
    optimizer <- keras3::optimizer_adam(learning_rate = learningRate)

    # ========== TRAINING STEP ==========
    trainStep <- function(xBatch, currentLambda) {
        with(tf$GradientTape() %as% tape, {
            # Encode
            z <- encoder(xBatch, training = TRUE)

            # Decode
            xReconstructed <- decoder(z, training = TRUE)

            # ========== LOSSES ==========
            # 1. Reconstruction loss (MSE)
            reconLoss <- tf$reduce_mean(tf$square(xBatch - xReconstructed)) *
                         tf$cast(inputDim, tf$float32)

            # 2. Sliced Wasserstein Distance to N(0,1) prior
            priorSamples <- tf$random$normal(shape = tf$shape(z))
            swLoss <- slicedWassersteinDistance(z, priorSamples, numProjections)

            # Total loss (using currentLambda for annealing)
            totalLoss <- reconLoss + currentLambda * swLoss
        })

        gradients <- tape$gradient(totalLoss, allVariables)
        optimizer$apply_gradients(Map(list, gradients, allVariables))

        list(total = totalLoss, recon = reconLoss, sw = swLoss)
    }

    # ========== TRAINING LOOP ==========
    message("\nStarting training...")

    nBatches <- ceiling(nTrain / batchSize)
    history <- list(total = c(), recon = c(), sw = c(), lambda = c())

    for (epoch in seq_len(epochs)) {
        # Lambda annealing: linearly increase from 0 to lambda over lambdaAnnealEpochs
        if (lambdaAnnealEpochs > 0 && epoch <= lambdaAnnealEpochs) {
            currentLambda <- lambda * (epoch / lambdaAnnealEpochs)
        } else {
            currentLambda <- lambda
        }

        # Shuffle
        idx <- sample(nTrain)
        xShuffled <- xTrain[idx, , drop = FALSE]

        epochLosses <- list(total = 0, recon = 0, sw = 0)

        for (batch in seq_len(nBatches)) {
            startIdx <- (batch - 1) * batchSize + 1
            endIdx <- min(batch * batchSize, nTrain)
            batchIdx <- startIdx:endIdx

            xBatch <- tf$constant(xShuffled[batchIdx, , drop = FALSE], dtype = tf$float32)

            losses <- trainStep(xBatch, currentLambda)

            epochLosses$total <- epochLosses$total + as.numeric(losses$total)
            epochLosses$recon <- epochLosses$recon + as.numeric(losses$recon)
            epochLosses$sw <- epochLosses$sw + as.numeric(losses$sw)
        }

        # Average losses
        epochLosses <- lapply(epochLosses, function(x) x / nBatches)

        history$total <- c(history$total, epochLosses$total)
        history$recon <- c(history$recon, epochLosses$recon)
        history$sw <- c(history$sw, epochLosses$sw)
        history$lambda <- c(history$lambda, currentLambda)

        if (verbose > 0 && (epoch %% verbose == 0 || epoch == 1)) {
            lambdaInfo <- if (lambdaAnnealEpochs > 0 && epoch <= lambdaAnnealEpochs) {
                sprintf(", \u03bb: %.2f", currentLambda)
            } else { "" }

            message(sprintf("Epoch %3d/%d - loss: %.4f (recon: %.4f, sw: %.4f%s)",
                          epoch, epochs, epochLosses$total, epochLosses$recon,
                          epochLosses$sw, lambdaInfo))
        }
    }

    # Compute validation loss if provided
    valLoss <- NULL
    if (!is.null(valData)) {
        xVal <- tf$constant(as.matrix(valData[, useMarkers]), dtype = tf$float32)
        zVal <- encoder(xVal, training = FALSE)
        xValRecon <- decoder(zVal, training = FALSE)
        valLoss <- as.numeric(tf$reduce_mean(tf$square(xVal - xValRecon)))
        message(sprintf("\nValidation reconstruction MSE: %.6f", valLoss))
    }

    message("\nTraining complete!")

    # ========== RETURN MODEL ==========
    list(
        encoder = encoder,
        decoder = decoder,
        latentDim = latentDim,
        inputDim = inputDim,
        history = history,
        valLoss = valLoss,
        type = "swae_enhanced"
    )
}


#' Encode Data with Enhanced SWAE
#'
#' Get latent representations from a trained enhanced SWAE model.
#'
#' @param model Trained enhanced SWAE model from trainSWAEEnhanced.
#' @param newData Matrix or data frame of data to encode.
#' @return Matrix of latent representations (cells x latent dimensions).
#' @export
encodeSWAE <- function(model, newData) {
    tf <- tensorflow::tf
    newData <- as.matrix(newData)
    xTensor <- tf$constant(newData, dtype = tf$float32)

    latent <- model$encoder(xTensor, training = FALSE)
    latentMatrix <- as.matrix(latent)
    colnames(latentMatrix) <- paste0("latent_", seq_len(ncol(latentMatrix)))

    latentMatrix
}


#' Decode Latent Representations with Enhanced SWAE
#'
#' Reconstruct data from latent representations using trained enhanced SWAE decoder.
#'
#' @param model Trained enhanced SWAE model from trainSWAEEnhanced.
#' @param latentData Matrix of latent representations.
#' @return Matrix of reconstructed data.
#' @export
decodeSWAE <- function(model, latentData) {
    tf <- tensorflow::tf
    latentData <- as.matrix(latentData)
    zTensor <- tf$constant(latentData, dtype = tf$float32)

    reconstructed <- model$decoder(zTensor, training = FALSE)
    as.matrix(reconstructed)
}


#' Train SWAE with Clustering Head
#'
#' Trains a Sliced Wasserstein Autoencoder with an additional clustering head
#' for joint batch correction and clustering. The model learns to reconstruct
#' data while simultaneously learning cluster assignments.
#'
#' @param trainData Data frame or matrix of training data.
#' @param useMarkers Character vector; names of marker columns to use.
#' @param nClusters Integer; number of clusters. Default is 10.
#' @param pseudoLabels Optional factor or integer vector of cluster labels for supervision.
#'   If provided, uses supervised cross-entropy loss. If NULL, uses DEC-style self-training.
#' @param epochs Integer; number of training epochs. Default is 100.
#' @param latentDim Integer; dimensionality of latent space. Default is 16.
#' @param lambda Numeric; weight for Sliced Wasserstein loss. Default is 0.1.
#' @param lambdaCluster Numeric; weight for clustering loss. Default is 1.0.
#' @param numProjections Integer; number of random projections for SW distance. Default is 50.
#' @param batchSize Integer; batch size. Default is 32.
#' @param learningRate Numeric; learning rate. Default is 1e-3.
#' @param hiddenSizes Integer vector; sizes of hidden layers. Default is c(23, 19).
#' @param clusterWarmup Integer; epochs before adding clustering loss. Default is 10.
#'   Allows the encoder to learn basic representations before clustering.
#' @param jointFineTune Logical; if TRUE, adds Stage 3 joint fine-tuning where
#'   encoder, decoder, and clustering head are all trained together. Default is FALSE.
#' @param jointEpochs Integer; number of epochs for joint fine-tuning. Default is 30.
#' @param labelSmoothing Numeric; label smoothing factor (0-1). Default is 0.0 (disabled).
#'   Values like 0.1 can reduce overconfidence and improve generalization.
#' @param contrastiveLoss Logical; if TRUE, adds supervised contrastive loss in Stage 3
#'   to push same-cluster samples together and different-cluster samples apart. Default is FALSE.
#' @param lambdaContrastive Numeric; weight for contrastive loss. Default is 0.1.
#' @param contrastiveTemp Numeric; temperature for contrastive loss (lower = sharper). Default is 0.5.
#' @param useMMD Logical; if TRUE, adds MMD loss to the training. Default is FALSE.
#' @param lambdaMMD Numeric; weight for MMD loss. Default is 0.1.
#' @param mmdBatchAlign Logical; if TRUE, adds MMD batch alignment loss. Default is FALSE.
#' @param lambdaMmdBatch Numeric; weight for MMD batch alignment loss. Default is 0.1.
#' @param batchCol Character; column name for batch IDs. Default is "sample_id".
#' @param contrastivePretraining Logical; if TRUE, adds contrastive pretraining stage. Default is FALSE.
#' @param contrastivePretrainEpochs Integer; number of epochs for contrastive pretraining. Default is 20.
#' @param lambdaContrastivePretrain Numeric; weight for contrastive pretraining loss. Default is 0.5.
#' @param augNoiseScale Numeric; noise scale for data augmentation. Default is 0.1.
#' @param augDropoutRate Numeric; dropout rate for data augmentation. Default is 0.1.
#' @param dropout Numeric; dropout rate for encoder/decoder layers (0-1). Default is 0.0 (disabled).
#' @param useBatchNorm Logical; if TRUE, add batch normalization layers. Default is FALSE.
#' @param useLayerNorm Logical; if TRUE, add layer normalization instead of batch norm. Default is FALSE.
#' @param weightDecay Numeric; L2 regularization weight for optimizer. Default is 0.0 (disabled).
#' @param lrSchedule Character; learning rate schedule type: "none", "cosine", "exponential". Default is "none".
#' @param lrDecaySteps Integer; steps for learning rate decay (for exponential). Default is 1000.
#' @param lrDecayRate Numeric; decay rate for exponential schedule. Default is 0.96.
#' @param earlyStoppingPatience Integer; epochs without improvement before stopping. 0 = disabled. Default is 0.
#' @param valSplit Numeric; fraction of training data to use for validation if valData not provided. Default is 0.1.
#' @param useContractive Logical; if TRUE, adds contractive regularization (DeepCAE). Default is FALSE.
#' @param lambdaContractive Numeric; weight for contractive loss. Default is 0.1.
#' @param useSCARF Logical; if TRUE, enables SCARF-style pretraining. Default is FALSE.
#' @param scarfCorruptionRate Numeric; fraction of features to corrupt in SCARF (0-1). Default is 0.3.
#' @param lambdaSCARF Numeric; weight for SCARF contrastive loss. Default is 1.0.
#' @param useAuxiliaryTask Logical; if TRUE, adds sample-level classification as auxiliary task. Default is FALSE.
#' @param lambdaAuxiliary Numeric; weight for auxiliary task loss. Default is 0.5.
#' @param auxiliaryLabels Vector of sample-level labels for auxiliary task.
#' @param auxiliarySampleCol Character; column name for sample IDs. Default is "sample_id".
#' @param useSwitchTab Logical; if TRUE, enables SwitchTab-style salient/mutual feature decoupling. Default is FALSE.
#' @param salientDim Integer; dimensionality of salient features. Default is NULL (same as latentDim).
#' @param lambdaSwitchTab Numeric; weight for SwitchTab losses. Default is 1.0.
#' @param valData Optional validation data.
#' @param seed Integer; random seed for reproducibility. Default is 1994.
#' @return List containing encoder, decoder, clusterHead, training history, and optionally
#'   salientEncoder (if useSwitchTab) and auxiliaryHead (if useAuxiliaryTask).
#' @importFrom tensorflow tf set_random_seed
#' @importFrom keras3 layer_input layer_dense keras_model
#' @export
#' @examples
#' \dontrun{
#' # With pseudo-labels from FuseSOM
#' model <- trainSWAEClustering(
#'   trainData = train_data,
#'   useMarkers = markers,
#'   nClusters = 11,
#'   pseudoLabels = sce$clusters,
#'   epochs = 100,
#'   latentDim = 16,
#'   lambda = 0.1,
#'   lambdaCluster = 1.0
#' )
#'
#' # Batch correction via decode
#' latent <- as.matrix(model$encoder(as.matrix(data[, markers])))
#' corrected <- as.matrix(model$decoder(latent))
#'
#' # Get cluster assignments
#' clusters <- as.matrix(model$clusterHead(latent))
#' cluster_labels <- max.col(clusters)
#' }
trainSWAEClustering <- function(trainData, useMarkers, nClusters = 10L,
                                 pseudoLabels = NULL, epochs = 100L, latentDim = 16L,
                                 lambda = 0.1, lambdaCluster = 1.0, numProjections = 50L,
                                 batchSize = 32L, learningRate = 1e-3,
                                 hiddenSizes = c(23L, 19L), clusterWarmup = 10L,
                                 jointFineTune = FALSE, jointEpochs = 30L,
                                 labelSmoothing = 0.0, contrastiveLoss = FALSE,
                                 lambdaContrastive = 0.1, contrastiveTemp = 0.5,
                                 useMMD = FALSE, lambdaMMD = 0.1,
                                 mmdBatchAlign = FALSE, lambdaMmdBatch = 0.1, batchCol = "sample_id",
                                 contrastivePretraining = FALSE, contrastivePretrainEpochs = 20L,
                                 lambdaContrastivePretrain = 0.5, augNoiseScale = 0.1,
                                 augDropoutRate = 0.1,
                                 dropout = 0.0, useBatchNorm = FALSE, useLayerNorm = FALSE,
                                 weightDecay = 0.0, lrSchedule = "none",
                                 lrDecaySteps = 1000L, lrDecayRate = 0.96,
                                 earlyStoppingPatience = 0L, valSplit = 0.1,
                                 # New representation learning methods
                                 useContractive = FALSE, lambdaContractive = 0.1,
                                 useSCARF = FALSE, scarfCorruptionRate = 0.3, lambdaSCARF = 1.0,
                                 useAuxiliaryTask = FALSE, lambdaAuxiliary = 0.5,
                                 auxiliaryLabels = NULL, auxiliarySampleCol = "sample_id",
                                 useSwitchTab = FALSE, salientDim = NULL, lambdaSwitchTab = 1.0,
                                 useDSW = FALSE, dswMaxIter = 10L, dswLambda = 1.0,
                                 valData = NULL, seed = 1994) {

    tensorflow::set_random_seed(seed)
    tf <- tensorflow::tf
    set.seed(seed)

    # Prepare data
    xAll <- as.matrix(trainData[, useMarkers])
    inputDim <- ncol(xAll)
    nAll <- nrow(xAll)

    # Handle validation split
    xVal <- NULL
    yVal <- NULL
    yValOneHot <- NULL
    useEarlyStopping <- earlyStoppingPatience > 0

    if (!is.null(valData)) {
        xVal <- as.matrix(valData[, useMarkers])
        xTrain <- xAll
        nTrain <- nAll
    } else if (useEarlyStopping && valSplit > 0) {
        # Split training data
        nVal <- round(nAll * valSplit)
        valIdx <- sample(nAll, nVal)
        trainIdx <- setdiff(seq_len(nAll), valIdx)
        xVal <- xAll[valIdx, , drop = FALSE]
        xTrain <- xAll[trainIdx, , drop = FALSE]
        nTrain <- nrow(xTrain)
        message(sprintf("  Validation split: %d train, %d val (%.0f%%)",
                       nTrain, nVal, valSplit * 100))
    } else {
        xTrain <- xAll
        nTrain <- nAll
    }

    # Handle batch alignment
    batchIndices <- NULL
    uniqueBatches <- NULL
    if (mmdBatchAlign) {
        if (!batchCol %in% colnames(trainData)) {
            warning(sprintf("Batch column '%s' not found, disabling batch alignment", batchCol))
            mmdBatchAlign <- FALSE
        } else {
            batchIds <- trainData[[batchCol]]
            uniqueBatches <- unique(batchIds)
            # Create list of indices for each batch
            batchIndices <- lapply(uniqueBatches, function(b) which(batchIds == b))
            names(batchIndices) <- uniqueBatches
            message(sprintf("  Batch alignment: %d batches found", length(uniqueBatches)))
        }
    }

    # Handle pseudo-labels
    supervised <- !is.null(pseudoLabels)
    if (supervised) {
        # Convert to factor first to handle character labels like "cluster_1"
        if (!is.factor(pseudoLabels)) {
            pseudoLabels <- as.factor(pseudoLabels)
        }
        # Get numeric codes (1-indexed from factor)
        yAll <- as.integer(pseudoLabels) - 1L  # Convert to 0-indexed

        # Check for valid labels
        if (any(is.na(yAll))) {
            stop("Invalid pseudo-labels: could not convert to integers")
        }
        if (max(yAll) >= nClusters) {
            stop(sprintf("Pseudo-labels have %d unique values but nClusters=%d",
                        max(yAll) + 1, nClusters))
        }

        # Split labels if validation split was done
        if (exists("valIdx") && !is.null(valIdx)) {
            yVal <- yAll[valIdx]
            yTrain <- yAll[trainIdx]
        } else {
            yTrain <- yAll
        }

        message(sprintf("  Label range: %d to %d (expecting 0 to %d)",
                       min(yTrain), max(yTrain), nClusters - 1))

        # One-hot encode training labels
        yTrainOneHot <- matrix(0, nrow = nTrain, ncol = nClusters)
        for (i in seq_len(nTrain)) {
            yTrainOneHot[i, yTrain[i] + 1] <- 1
        }

        # One-hot encode validation labels if present
        if (!is.null(yVal)) {
            yValOneHot <- matrix(0, nrow = length(yVal), ncol = nClusters)
            for (i in seq_len(length(yVal))) {
                yValOneHot[i, yVal[i] + 1] <- 1
            }
        }

        message(sprintf("  One-hot matrix: %d x %d, sum per row: %.1f",
                       nrow(yTrainOneHot), ncol(yTrainOneHot), mean(rowSums(yTrainOneHot))))
    }

    # Ensure integer types
    nClusters <- as.integer(nClusters)
    latentDim <- as.integer(latentDim)
    batchSize <- as.integer(batchSize)
    numProjections <- as.integer(numProjections)
    hiddenSizes <- as.integer(hiddenSizes)
    clusterWarmup <- as.integer(clusterWarmup)

    intermediateDim <- hiddenSizes[1]
    intermediateDim2 <- hiddenSizes[2]

    message("SWAE-Clustering Architecture:")
    message(sprintf("  Input dim: %d", inputDim))
    message(sprintf("  Hidden: %d -> %d -> latent(%d)", intermediateDim, intermediateDim2, latentDim))
    message(sprintf("  Clusters: %d", nClusters))
    message(sprintf("  Training: %s mode", ifelse(supervised, "supervised", "self-supervised")))
    message(sprintf("  Lambda (SW): %.2f, Lambda (cluster): %.2f", lambda, lambdaCluster))
    message(sprintf("  Cluster warmup: %d epochs", clusterWarmup))
    if (useMMD) {
        message(sprintf("  MMD loss: enabled (lambda=%.2f)", lambdaMMD))
    }
    if (mmdBatchAlign) {
        message(sprintf("  MMD batch alignment: enabled (lambda=%.2f, %d batches)",
                       lambdaMmdBatch, length(uniqueBatches)))
    }
    if (dropout > 0) {
        message(sprintf("  Dropout: %.2f", dropout))
    }
    if (useBatchNorm) {
        message("  Batch normalization: enabled")
    }
    if (useLayerNorm) {
        message("  Layer normalization: enabled")
    }
    if (weightDecay > 0) {
        message(sprintf("  Weight decay: %.6f", weightDecay))
    }
    if (lrSchedule != "none") {
        message(sprintf("  LR schedule: %s", lrSchedule))
    }
    if (useEarlyStopping) {
        message(sprintf("  Early stopping: patience=%d", earlyStoppingPatience))
    }
    if (useContractive) {
        message(sprintf("  Contractive regularization (DeepCAE): lambda=%.4f", lambdaContractive))
    }
    if (useSCARF) {
        message(sprintf("  SCARF pretraining: corruption=%.2f, lambda=%.2f", scarfCorruptionRate, lambdaSCARF))
    }
    if (useAuxiliaryTask) {
        message(sprintf("  Auxiliary task supervision: lambda=%.2f", lambdaAuxiliary))
    }
    if (useSwitchTab) {
        salientDimActual <- if (is.null(salientDim)) latentDim else salientDim
        message(sprintf("  SwitchTab (salient/mutual): salient_dim=%d, lambda=%.2f", salientDimActual, lambdaSwitchTab))
    }
    if (useDSW) {
        message(sprintf("  Distributional SW: max_iter=%d, lambda=%.2f", dswMaxIter, dswLambda))
    }

    # ========== BUILD ENCODER ==========
    encoderInput <- keras3::layer_input(shape = inputDim, name = "encoder_input")
    h <- encoderInput |>
        keras3::layer_dense(intermediateDim, activation = "relu", name = "enc_h1")
    if (useBatchNorm) {
        h <- h |> keras3::layer_batch_normalization(name = "enc_bn1")
    }
    if (useLayerNorm) {
        h <- h |> keras3::layer_layer_normalization(name = "enc_ln1")
    }
    if (dropout > 0) {
        h <- h |> keras3::layer_dropout(rate = dropout, name = "enc_drop1")
    }
    h <- h |>
        keras3::layer_dense(intermediateDim2, activation = "relu", name = "enc_h2")
    if (useBatchNorm) {
        h <- h |> keras3::layer_batch_normalization(name = "enc_bn2")
    }
    if (useLayerNorm) {
        h <- h |> keras3::layer_layer_normalization(name = "enc_ln2")
    }
    if (dropout > 0) {
        h <- h |> keras3::layer_dropout(rate = dropout, name = "enc_drop2")
    }
    h <- h |>
        keras3::layer_dense(latentDim, activation = "gelu", name = "latent")
    encoder <- keras3::keras_model(encoderInput, h, name = "encoder")

    # ========== BUILD DECODER ==========
    decoderInput <- keras3::layer_input(shape = latentDim, name = "decoder_input")
    decoded <- decoderInput |>
        keras3::layer_dense(intermediateDim2, activation = "relu", name = "dec_h1")
    if (useBatchNorm) {
        decoded <- decoded |> keras3::layer_batch_normalization(name = "dec_bn1")
    }
    if (useLayerNorm) {
        decoded <- decoded |> keras3::layer_layer_normalization(name = "dec_ln1")
    }
    if (dropout > 0) {
        decoded <- decoded |> keras3::layer_dropout(rate = dropout, name = "dec_drop1")
    }
    decoded <- decoded |>
        keras3::layer_dense(intermediateDim, activation = "relu", name = "dec_h2")
    if (useBatchNorm) {
        decoded <- decoded |> keras3::layer_batch_normalization(name = "dec_bn2")
    }
    if (useLayerNorm) {
        decoded <- decoded |> keras3::layer_layer_normalization(name = "dec_ln2")
    }
    if (dropout > 0) {
        decoded <- decoded |> keras3::layer_dropout(rate = dropout, name = "dec_drop2")
    }
    decoded <- decoded |>
        keras3::layer_dense(inputDim, activation = "sigmoid", name = "reconstruction")
    decoder <- keras3::keras_model(decoderInput, decoded, name = "decoder")

    # ========== BUILD CLUSTERING HEAD ==========
    clusterInput <- keras3::layer_input(shape = latentDim, name = "cluster_input")
    clusterLogits <- clusterInput |>
        keras3::layer_dense(32L, activation = "relu", name = "cluster_h1")
    if (dropout > 0) {
        clusterLogits <- clusterLogits |> keras3::layer_dropout(rate = dropout / 2, name = "cluster_drop")
    }
    clusterLogits <- clusterLogits |>
        keras3::layer_dense(nClusters, activation = "softmax", name = "cluster_output")
    clusterHead <- keras3::keras_model(clusterInput, clusterLogits, name = "cluster_head")

    # ========== BUILD SWITCHTAB SALIENT ENCODER (optional) ==========
    salientEncoder <- NULL
    if (useSwitchTab) {
        salientDimActual <- if (is.null(salientDim)) latentDim else as.integer(salientDim)
        # Salient encoder: maps latent to salient features for downstream tasks
        salientInput <- keras3::layer_input(shape = latentDim, name = "salient_input")
        salientH <- salientInput |>
            keras3::layer_dense(latentDim, activation = "relu", name = "salient_h1") |>
            keras3::layer_dense(salientDimActual, activation = "linear", name = "salient_output")
        salientEncoder <- keras3::keras_model(salientInput, salientH, name = "salient_encoder")

        # Mutual encoder: maps latent to mutual features for reconstruction
        mutualInput <- keras3::layer_input(shape = latentDim, name = "mutual_input")
        mutualH <- mutualInput |>
            keras3::layer_dense(latentDim, activation = "relu", name = "mutual_h1") |>
            keras3::layer_dense(latentDim - salientDimActual, activation = "linear", name = "mutual_output")
        mutualEncoder <- keras3::keras_model(mutualInput, mutualH, name = "mutual_encoder")

        # Combined decoder for SwitchTab (takes concatenated salient + mutual)
        switchDecoderInput <- keras3::layer_input(shape = latentDim, name = "switch_decoder_input")
        switchDecoded <- switchDecoderInput |>
            keras3::layer_dense(intermediateDim2, activation = "relu", name = "switch_dec_h1") |>
            keras3::layer_dense(intermediateDim, activation = "relu", name = "switch_dec_h2") |>
            keras3::layer_dense(inputDim, activation = "sigmoid", name = "switch_reconstruction")
        switchDecoder <- keras3::keras_model(switchDecoderInput, switchDecoded, name = "switch_decoder")
    }

    # ========== BUILD AUXILIARY TASK HEAD (optional) ==========
    auxiliaryHead <- NULL
    if (useAuxiliaryTask) {
        # Prepare sample-level labels
        if (is.null(auxiliaryLabels)) {
            stop("auxiliaryLabels required when useAuxiliaryTask=TRUE")
        }
        # Get unique sample labels
        sampleIds <- trainData[[auxiliarySampleCol]]
        uniqueSamples <- unique(sampleIds)
        sampleLabelMap <- setNames(auxiliaryLabels[match(uniqueSamples, sampleIds)], uniqueSamples)

        # Build auxiliary prediction head (binary classification)
        auxInput <- keras3::layer_input(shape = latentDim, name = "aux_input")
        auxH <- auxInput |>
            keras3::layer_dense(32L, activation = "relu", name = "aux_h1") |>
            keras3::layer_dense(1L, activation = "sigmoid", name = "aux_output")
        auxiliaryHead <- keras3::keras_model(auxInput, auxH, name = "auxiliary_head")
    }

    # ========== BUILD DSW PROJECTION NETWORK (optional) ==========
    dswProjectionNet <- NULL
    dswProjectionOptimizer <- NULL
    if (useDSW) {
        # Small network that transforms random projections to optimized ones
        # Input: [numProjections, latentDim] -> Output: [numProjections, latentDim]
        dswInput <- keras3::layer_input(shape = latentDim, name = "dsw_input")
        dswH <- dswInput |>
            keras3::layer_dense(latentDim, activation = "relu", name = "dsw_h1") |>
            keras3::layer_dense(latentDim, activation = "linear", name = "dsw_output")
        dswProjectionNet <- keras3::keras_model(dswInput, dswH, name = "dsw_projection")
        dswProjectionOptimizer <- keras3::optimizer_adam(learning_rate = 1e-4)
        message("  DSW projection network built")
    }

    # ========== CONTRASTIVE LOSS FUNCTION ==========
    # Supervised contrastive loss: pulls same-class samples together, pushes different apart
    supervisedContrastiveLoss <- function(z, labels, temperature = 0.5) {
        # Normalize embeddings
        zNorm <- tf$math$l2_normalize(z, axis = 1L)

        # Compute similarity matrix
        similarity <- tf$matmul(zNorm, zNorm, transpose_b = TRUE) / temperature

        # Create mask for positive pairs (same label)
        labelsExpanded <- tf$expand_dims(labels, 1L)
        positiveMask <- tf$cast(tf$equal(labelsExpanded, tf$transpose(labelsExpanded)), tf$float32)

        # Remove diagonal (self-similarity)
        batchSize <- tf$shape(z)[1]
        diagMask <- tf$ones(c(batchSize, batchSize)) - tf$eye(batchSize)
        positiveMask <- positiveMask * diagMask

        # Compute log softmax
        logits <- similarity - tf$reduce_max(similarity, axis = 1L, keepdims = TRUE)
        expLogits <- tf$exp(logits) * diagMask
        logSumExp <- tf$math$log(tf$reduce_sum(expLogits, axis = 1L, keepdims = TRUE) + 1e-10)
        logProb <- logits - logSumExp

        # Average over positive pairs
        nPositives <- tf$reduce_sum(positiveMask, axis = 1L)
        nPositives <- tf$maximum(nPositives, 1.0)  # Avoid division by zero
        meanLogProb <- tf$reduce_sum(positiveMask * logProb, axis = 1L) / nPositives

        loss <- -tf$reduce_mean(meanLogProb)
        return(loss)
    }

    # ========== CONTRACTIVE LOSS (DeepCAE) ==========
    computeContractiveLoss <- function(xBatch) {
        xBatchTensor <- tf$constant(xBatch, dtype = tf$float32)
        with(tf$GradientTape(persistent = TRUE) %as% tape, {
            tape$watch(xBatchTensor)
            z <- encoder(xBatchTensor, training = TRUE)
        })
        jacobianNorm <- tf$constant(0.0)
        for (i in seq_len(latentDim)) {
            grad <- tape$gradient(z[, i], xBatchTensor)
            if (!is.null(grad)) {
                jacobianNorm <- jacobianNorm + tf$reduce_mean(tf$reduce_sum(tf$square(grad), axis = 1L))
            }
        }
        return(jacobianNorm / tf$cast(latentDim, tf$float32))
    }

    # ========== SCARF LOSS (Random Feature Corruption) ==========
    scarfCorruptFeatures <- function(xBatch) {
        batchSizeTF <- tf$shape(xBatch)[1]
        nFeatures <- tf$shape(xBatch)[2]

        mask <- tf$cast(tf$random$uniform(tf$shape(xBatch)) < scarfCorruptionRate, tf$float32)
        randomIdx <- tf$random$shuffle(tf$range(batchSizeTF))
        xShuffled <- tf$gather(xBatch, randomIdx)
        xCorrupted <- xBatch * (1.0 - mask) + xShuffled * mask
        return(xCorrupted)
    }

    scarfContrastiveLoss <- function(z1, z2, temperature = 0.5) {
        batchSizeTF <- tf$shape(z1)[1]
        z1Norm <- tf$math$l2_normalize(z1, axis = 1L)
        z2Norm <- tf$math$l2_normalize(z2, axis = 1L)
        posSim <- tf$reduce_sum(z1Norm * z2Norm, axis = 1L) / temperature
        allSim <- tf$matmul(z1Norm, z2Norm, transpose_b = TRUE) / temperature
        labels <- tf$range(batchSizeTF)
        loss <- tf$reduce_mean(tf$nn$sparse_softmax_cross_entropy_with_logits(
            labels = labels, logits = allSim
        ))
        return(loss)
    }

    # ========== SWITCHTAB LOSS ==========
    switchTabLoss <- function(z1, z2, xBatch1, xBatch2) {
        salient1 <- salientEncoder(z1, training = TRUE)
        salient2 <- salientEncoder(z2, training = TRUE)
        mutual1 <- mutualEncoder(z1, training = TRUE)
        mutual2 <- mutualEncoder(z2, training = TRUE)

        combined12 <- tf$concat(list(salient2, mutual1), axis = 1L)
        combined21 <- tf$concat(list(salient1, mutual2), axis = 1L)
        recon12 <- switchDecoder(combined12, training = TRUE)
        recon21 <- switchDecoder(combined21, training = TRUE)

        combined11 <- tf$concat(list(salient1, mutual1), axis = 1L)
        combined22 <- tf$concat(list(salient2, mutual2), axis = 1L)
        recon11 <- switchDecoder(combined11, training = TRUE)
        recon22 <- switchDecoder(combined22, training = TRUE)

        selfReconLoss <- (tf$reduce_mean(tf$square(xBatch1 - recon11)) +
                          tf$reduce_mean(tf$square(xBatch2 - recon22))) / 2.0
        salientContrLoss <- scarfContrastiveLoss(salient1, salient2, temperature = 0.5)
        return(list(selfRecon = selfReconLoss, salientContr = salientContrLoss))
    }

    if (contrastiveLoss && supervised) {
        message(sprintf("  Contrastive loss: enabled (lambda=%.2f, temp=%.2f)",
                       lambdaContrastive, contrastiveTemp))
    }
    if (labelSmoothing > 0) {
        message(sprintf("  Label smoothing: %.2f", labelSmoothing))
    }
    if (contrastivePretraining) {
        message(sprintf("  Contrastive pretraining: %d epochs (lambda=%.2f, noise=%.2f, dropout=%.2f)",
                       contrastivePretrainEpochs, lambdaContrastivePretrain, augNoiseScale, augDropoutRate))
    }

    # ========== STAGE 0: CONTRASTIVE PRETRAINING (optional) ==========
    if (contrastivePretraining) {
        message("\n=== STAGE 0: Contrastive pretraining (SimCLR-style) ===")

        optimizerPretrain <- keras3::optimizer_adam(learning_rate = learningRate)
        encoderVars <- encoder$trainable_variables

        contrastivePretrainStep <- function(xBatch) {
            with(tf$GradientTape() %as% tape, {
                xAug1 <- augmentCells(xBatch, noiseScale = augNoiseScale,
                                      dropoutRate = augDropoutRate)
                xAug2 <- augmentCells(xBatch, noiseScale = augNoiseScale,
                                      dropoutRate = augDropoutRate)
                z1 <- encoder(xAug1, training = TRUE)
                z2 <- encoder(xAug2, training = TRUE)
                contrLoss <- unsupervisedContrastiveLoss(z1, z2, temperature = contrastiveTemp)
                zOrig <- encoder(xBatch, training = TRUE)
                xRecon <- decoder(zOrig, training = TRUE)
                reconLoss <- tf$reduce_mean(tf$square(xBatch - xRecon)) *
                    tf$cast(inputDim, tf$float32)
                totalLoss <- reconLoss + lambdaContrastivePretrain * contrLoss
            })
            allVars <- c(encoder$trainable_variables, decoder$trainable_variables)
            gradients <- tape$gradient(totalLoss, allVars)
            optimizerPretrain$apply_gradients(Map(list, gradients, allVars))
            list(total = totalLoss, contrastive = contrLoss, recon = reconLoss)
        }

        nBatchesPretrain <- ceiling(nTrain / batchSize)
        for (epoch in seq_len(contrastivePretrainEpochs)) {
            idx <- sample(nTrain)
            xShuffled <- xTrain[idx, , drop = FALSE]

            epochLoss <- 0
            epochContr <- 0
            for (batch in seq_len(nBatchesPretrain)) {
                startIdx <- (batch - 1) * batchSize + 1
                endIdx <- min(batch * batchSize, nTrain)
                xBatch <- xShuffled[startIdx:endIdx, , drop = FALSE]
                losses <- contrastivePretrainStep(xBatch)
                epochLoss <- epochLoss + as.numeric(losses$total)
                epochContr <- epochContr + as.numeric(losses$contrastive)
            }

            if (epoch %% 5 == 0 || epoch == 1) {
                message(sprintf("Stage 0 - Epoch %d/%d - total: %.4f, contrastive: %.4f",
                               epoch, contrastivePretrainEpochs,
                               epochLoss / nBatchesPretrain,
                               epochContr / nBatchesPretrain))
            }
        }
        message("Contrastive pretraining complete!")
    }

    # ========== STAGE 1: TRAIN SWAE ==========
    if (useMMD) {
        message("\n=== STAGE 1: Training SWAE (reconstruction + SW + MMD) ===")
    } else {
        message("\n=== STAGE 1: Training SWAE (reconstruction + SW) ===")
    }

    # Create learning rate schedule
    nBatchesTotal <- ceiling(nTrain / batchSize)
    totalSteps <- nBatchesTotal * epochs
    if (lrSchedule == "cosine") {
        lrScheduleObj <- tf$keras$optimizers$schedules$CosineDecay(
            initial_learning_rate = learningRate,
            decay_steps = as.integer(totalSteps),
            alpha = 0.01
        )
        currentLR <- lrScheduleObj
    } else if (lrSchedule == "exponential") {
        lrScheduleObj <- tf$keras$optimizers$schedules$ExponentialDecay(
            initial_learning_rate = learningRate,
            decay_steps = as.integer(lrDecaySteps),
            decay_rate = lrDecayRate
        )
        currentLR <- lrScheduleObj
    } else {
        currentLR <- learningRate
    }

    if (weightDecay > 0) {
        optimizerAE <- keras3::optimizer_adam_w(learning_rate = currentLR, weight_decay = weightDecay)
    } else {
        optimizerAE <- keras3::optimizer_adam(learning_rate = currentLR)
    }
    aeVariables <- c(encoder$trainable_variables, decoder$trainable_variables)

    # Early stopping state
    bestValLoss <- Inf
    patienceCounter <- 0
    bestEncoderWeights <- NULL
    bestDecoderWeights <- NULL

    swaeTrainStep <- function(xBatch, xBatch2 = NULL, sampleLabels = NULL) {
        xBatchArray <- as.array(as.matrix(xBatch))
        xBatchTensor <- tf$constant(xBatchArray, dtype = tf$float32)
        xBatch2Tensor <- NULL
        if (!is.null(xBatch2)) {
            xBatch2Array <- as.array(as.matrix(xBatch2))
            xBatch2Tensor <- tf$constant(xBatch2Array, dtype = tf$float32)
        }

        with(tf$GradientTape(persistent = TRUE) %as% tape, {
            tape$watch(xBatchTensor)
            z <- encoder(xBatchTensor, training = TRUE)
            xRecon <- decoder(z, training = TRUE)

            reconLoss <- tf$reduce_mean(tf$square(xBatchTensor - xRecon)) *
                tf$cast(inputDim, tf$float32)

            priorSamples <- tf$random$normal(shape = tf$shape(z))
            if (useDSW) {
                swLoss <- distributionalSlicedWassersteinDistance(
                    z, priorSamples, numProjections,
                    dswProjectionNet, dswProjectionOptimizer,
                    maxIter = dswMaxIter, lam = dswLambda)
            } else {
                swLoss <- slicedWassersteinDistance(z, priorSamples, numProjections)
            }

            totalLoss <- reconLoss + lambda * swLoss

            mmdLoss <- tf$constant(0.0)
            if (useMMD) {
                mmdLoss <- computeMMD(z, priorSamples)
                totalLoss <- totalLoss + lambdaMMD * mmdLoss
            }

            contractiveLossVal <- tf$constant(0.0)
            if (useContractive) {
                jacobianNorm <- tf$constant(0.0)
                for (i in seq_len(latentDim)) {
                    grad <- tape$gradient(z[, i], xBatchTensor)
                    if (!is.null(grad)) {
                        jacobianNorm <- jacobianNorm + tf$reduce_mean(tf$reduce_sum(tf$square(grad), axis = 1L))
                    }
                }
                contractiveLossVal <- jacobianNorm / tf$cast(latentDim, tf$float32)
                totalLoss <- totalLoss + lambdaContractive * contractiveLossVal
            }

            scarfLossVal <- tf$constant(0.0)
            if (useSCARF) {
                xCorrupted <- scarfCorruptFeatures(xBatchTensor)
                zCorrupted <- encoder(xCorrupted, training = TRUE)
                scarfLossVal <- scarfContrastiveLoss(z, zCorrupted, temperature = contrastiveTemp)
                totalLoss <- totalLoss + lambdaSCARF * scarfLossVal
            }

            switchLossVal <- tf$constant(0.0)
            if (useSwitchTab && !is.null(xBatch2Tensor)) {
                z2 <- encoder(xBatch2Tensor, training = TRUE)
                switchLosses <- switchTabLoss(z, z2, xBatchTensor, xBatch2Tensor)
                switchLossVal <- switchLosses$selfRecon + switchLosses$salientContr
                totalLoss <- totalLoss + lambdaSwitchTab * switchLossVal
            }

            auxLossVal <- tf$constant(0.0)
            if (useAuxiliaryTask && !is.null(sampleLabels)) {
                zMean <- tf$reduce_mean(z, axis = 0L, keepdims = TRUE)
                auxPred <- auxiliaryHead(zMean, training = TRUE)
                auxLabelMean <- tf$reduce_mean(tf$cast(sampleLabels, tf$float32))
                auxLossVal <- tf$keras$losses$binary_crossentropy(
                    tf$reshape(auxLabelMean, c(1L, 1L)),
                    auxPred
                )
                totalLoss <- totalLoss + lambdaAuxiliary * auxLossVal
            }
        })

        allVars <- aeVariables
        if (useSwitchTab && !is.null(salientEncoder)) {
            allVars <- c(allVars, salientEncoder$trainable_variables,
                        mutualEncoder$trainable_variables, switchDecoder$trainable_variables)
        }
        if (useAuxiliaryTask && !is.null(auxiliaryHead)) {
            allVars <- c(allVars, auxiliaryHead$trainable_variables)
        }

        gradients <- tape$gradient(totalLoss, allVars)
        optimizerAE$apply_gradients(Map(list, gradients, allVars))
        list(total = totalLoss, recon = reconLoss, sw = swLoss, mmd = mmdLoss,
             contractive = contractiveLossVal, scarf = scarfLossVal,
             switch = switchLossVal, auxiliary = auxLossVal)
    }

    # Batch alignment step
    batchAlignStep <- function(nSamplesPerBatch = 256L) {
        if (!mmdBatchAlign || length(uniqueBatches) < 2) return(tf$constant(0.0))

        batchPairIdx <- sample(length(uniqueBatches), 2)
        idx1 <- batchIndices[[batchPairIdx[1]]]
        idx2 <- batchIndices[[batchPairIdx[2]]]

        n1 <- min(nSamplesPerBatch, length(idx1))
        n2 <- min(nSamplesPerBatch, length(idx2))
        samp1 <- sample(idx1, n1)
        samp2 <- sample(idx2, n2)

        x1 <- xTrain[samp1, , drop = FALSE]
        x2 <- xTrain[samp2, , drop = FALSE]

        with(tf$GradientTape() %as% tape, {
            z1 <- encoder(x1, training = TRUE)
            z2 <- encoder(x2, training = TRUE)
            batchMmdLoss <- computeMMD(z1, z2)
        })

        gradients <- tape$gradient(batchMmdLoss, encoder$trainable_variables)
        gradients <- lapply(gradients, function(g) {
            if (!is.null(g)) tf$clip_by_norm(g, 1.0) else g
        })
        optimizerAE$apply_gradients(Map(list, gradients, encoder$trainable_variables))
        batchMmdLoss
    }

    # Validation loss computation function
    computeValLoss <- function() {
        if (is.null(xVal)) return(Inf)
        z <- encoder(xVal, training = FALSE)
        xRecon <- decoder(z, training = FALSE)
        reconLoss <- tf$reduce_mean(tf$square(xVal - xRecon)) * tf$cast(inputDim, tf$float32)
        as.numeric(reconLoss)
    }

    # Stage 1 training loop
    nBatches <- ceiling(nTrain / batchSize)
    stage1Epochs <- clusterWarmup

    # Prepare sample labels for auxiliary task if needed
    trainSampleIds <- NULL
    trainSampleLabels <- NULL
    if (useAuxiliaryTask && !is.null(auxiliaryLabels)) {
        trainSampleIds <- trainData[[auxiliarySampleCol]]
        uniqueSamples <- unique(trainSampleIds)
        sampleToLabel <- sapply(uniqueSamples, function(s) {
            auxiliaryLabels[which(trainSampleIds == s)[1]]
        })
        trainSampleLabels <- sampleToLabel[as.character(trainSampleIds)]
    }

    for (epoch in seq_len(stage1Epochs)) {
        idx <- sample(nTrain)
        xShuffled <- xTrain[idx, , drop = FALSE]

        labelsShuffled <- NULL
        if (useAuxiliaryTask && !is.null(trainSampleLabels)) {
            labelsShuffled <- trainSampleLabels[idx]
        }

        epochLoss <- 0
        epochBatchMmd <- 0
        epochContractive <- 0
        epochScarf <- 0
        epochSwitch <- 0
        epochAux <- 0

        for (batch in seq_len(nBatches)) {
            startIdx <- (batch - 1) * batchSize + 1
            endIdx <- min(batch * batchSize, nTrain)
            xBatch <- xShuffled[startIdx:endIdx, , drop = FALSE]

            xBatch2 <- NULL
            if (useSwitchTab) {
                idx2 <- sample(nTrain, endIdx - startIdx + 1)
                xBatch2 <- xTrain[idx2, , drop = FALSE]
            }

            batchLabels <- NULL
            if (useAuxiliaryTask && !is.null(labelsShuffled)) {
                batchLabels <- labelsShuffled[startIdx:endIdx]
            }

            losses <- swaeTrainStep(xBatch, xBatch2, batchLabels)
            epochLoss <- epochLoss + as.numeric(losses$total)
            if (useContractive) epochContractive <- epochContractive + as.numeric(losses$contractive)
            if (useSCARF) epochScarf <- epochScarf + as.numeric(losses$scarf)
            if (useSwitchTab) epochSwitch <- epochSwitch + as.numeric(losses$switch)
            if (useAuxiliaryTask) epochAux <- epochAux + as.numeric(losses$auxiliary)

            if (mmdBatchAlign && batch %% 5 == 0) {
                batchMmd <- batchAlignStep(256L)
                epochBatchMmd <- epochBatchMmd + as.numeric(batchMmd)
            }
        }

        valLoss <- computeValLoss()

        if (epoch %% 10 == 0 || epoch == 1) {
            logMsg <- sprintf("Stage 1 - Epoch %d/%d - loss: %.4f", epoch, stage1Epochs, epochLoss / nBatches)
            if (useEarlyStopping) {
                logMsg <- paste0(logMsg, sprintf(", val: %.4f", valLoss))
            }
            if (useContractive) {
                logMsg <- paste0(logMsg, sprintf(", contr: %.4f", epochContractive / nBatches))
            }
            if (useSCARF) {
                logMsg <- paste0(logMsg, sprintf(", scarf: %.4f", epochScarf / nBatches))
            }
            if (useSwitchTab) {
                logMsg <- paste0(logMsg, sprintf(", switch: %.4f", epochSwitch / nBatches))
            }
            if (useAuxiliaryTask) {
                logMsg <- paste0(logMsg, sprintf(", aux: %.4f", epochAux / nBatches))
            }
            if (mmdBatchAlign) {
                logMsg <- paste0(logMsg, sprintf(", batch_mmd: %.4f", epochBatchMmd / (nBatches / 5)))
            }
            message(logMsg)
        }

        if (useEarlyStopping && !is.null(xVal)) {
            if (valLoss < bestValLoss) {
                bestValLoss <- valLoss
                patienceCounter <- 0
                bestEncoderWeights <- encoder$get_weights()
                bestDecoderWeights <- decoder$get_weights()
            } else {
                patienceCounter <- patienceCounter + 1
                if (patienceCounter >= earlyStoppingPatience) {
                    message(sprintf("Early stopping at epoch %d (patience=%d)", epoch, earlyStoppingPatience))
                    if (!is.null(bestEncoderWeights)) {
                        encoder$set_weights(bestEncoderWeights)
                        decoder$set_weights(bestDecoderWeights)
                    }
                    break
                }
            }
        }
    }

    # ========== STAGE 2: TRAIN CLUSTERING HEAD ==========
    message("\n=== STAGE 2: Training clustering head (encoder frozen) ===")

    optimizerCluster <- keras3::optimizer_adam(learning_rate = learningRate)
    clusterVariables <- clusterHead$trainable_variables

    clusterTrainStep <- function(xBatch, yBatchOneHot) {
        z <- tf$stop_gradient(encoder(xBatch, training = FALSE))

        with(tf$GradientTape() %as% tape, {
            clusterProbs <- clusterHead(z, training = TRUE)

            if (supervised) {
                smoothedLabels <- yBatchOneHot * (1.0 - labelSmoothing) +
                    labelSmoothing / tf$cast(nClusters, tf$float32)
                clusterLoss <- tf$reduce_mean(
                    tf$keras$losses$categorical_crossentropy(smoothedLabels, clusterProbs)
                )
            } else {
                q <- clusterProbs + 1e-10
                f <- tf$reduce_sum(q, axis = 0L, keepdims = TRUE)
                p <- tf$square(q) / f
                p <- p / tf$reduce_sum(p, axis = 1L, keepdims = TRUE)
                clusterLoss <- tf$reduce_mean(
                    tf$reduce_sum(p * tf$math$log(p / q), axis = 1L)
                )
            }
        })

        gradients <- tape$gradient(clusterLoss, clusterVariables)
        optimizerCluster$apply_gradients(Map(list, gradients, clusterVariables))
        as.array(clusterLoss)
    }

    # Stage 2 training loop
    stage2Epochs <- epochs - clusterWarmup
    history <- list(cluster = c())

    for (epoch in seq_len(stage2Epochs)) {
        idx <- sample(nTrain)
        xShuffled <- xTrain[idx, , drop = FALSE]
        if (supervised) {
            yShuffled <- yTrainOneHot[idx, , drop = FALSE]
        }

        epochLoss <- 0
        for (batch in seq_len(nBatches)) {
            startIdx <- (batch - 1) * batchSize + 1
            endIdx <- min(batch * batchSize, nTrain)
            xBatch <- xShuffled[startIdx:endIdx, , drop = FALSE]

            if (supervised) {
                yBatch <- yShuffled[startIdx:endIdx, , drop = FALSE]
            } else {
                yBatch <- NULL
            }

            loss <- clusterTrainStep(xBatch, yBatch)
            epochLoss <- epochLoss + as.numeric(loss)
        }

        avgLoss <- epochLoss / nBatches
        history$cluster <- c(history$cluster, avgLoss)

        if (epoch %% 10 == 0 || epoch == 1) {
            message(sprintf("Stage 2 - Epoch %d/%d - cluster_loss: %.4f", epoch, stage2Epochs, avgLoss))
        }
    }

    # ========== STAGE 3: JOINT FINE-TUNING ==========
    if (supervised && (jointFineTune || contrastiveLoss)) {
        stageName <- if (contrastiveLoss && jointFineTune) {
            "Joint fine-tuning with contrastive loss"
        } else if (contrastiveLoss) {
            "Contrastive fine-tuning"
        } else {
            "Joint fine-tuning"
        }
        message(sprintf("\n=== STAGE 3: %s (all parameters) ===", stageName))

        optimizerJoint <- keras3::optimizer_adam(learning_rate = learningRate * 0.1, clipnorm = 1.0)
        allVariables <- c(encoder$trainable_variables, decoder$trainable_variables, clusterHead$trainable_variables)

        jointTrainStep <- function(xBatch, yBatchOneHot, yBatchLabels) {
            with(tf$GradientTape() %as% tape, {
                z <- encoder(xBatch, training = TRUE)
                xRecon <- decoder(z, training = TRUE)
                clusterProbs <- clusterHead(z, training = TRUE)

                reconLoss <- tf$reduce_mean(tf$square(xBatch - xRecon)) *
                    tf$cast(inputDim, tf$float32)

                priorSamples <- tf$random$normal(shape = tf$shape(z))
                if (useDSW) {
                    swLoss <- distributionalSlicedWassersteinDistance(
                        z, priorSamples, numProjections,
                        dswProjectionNet, dswProjectionOptimizer,
                        maxIter = dswMaxIter, lam = dswLambda)
                } else {
                    swLoss <- slicedWassersteinDistance(z, priorSamples, numProjections)
                }

                clusterLoss <- tf$reduce_mean(
                    tf$keras$losses$categorical_crossentropy(yBatchOneHot, clusterProbs)
                )

                totalLoss <- reconLoss + lambda * swLoss + lambdaCluster * clusterLoss

                mmdLoss <- tf$constant(0.0)
                if (useMMD) {
                    mmdLoss <- computeMMD(z, priorSamples)
                    totalLoss <- totalLoss + lambdaMMD * mmdLoss
                }

                contrLoss <- tf$constant(0.0)
                if (contrastiveLoss) {
                    contrLoss <- supervisedContrastiveLoss(z, yBatchLabels, contrastiveTemp)
                    totalLoss <- totalLoss + lambdaContrastive * contrLoss
                }
            })

            gradients <- tape$gradient(totalLoss, allVariables)
            gradients <- lapply(gradients, function(g) {
                if (!is.null(g)) tf$clip_by_norm(g, 1.0) else g
            })
            optimizerJoint$apply_gradients(Map(list, gradients, allVariables))
            list(total = as.array(totalLoss), recon = as.array(reconLoss),
                 sw = as.array(swLoss), cluster = as.array(clusterLoss),
                 contrastive = as.array(contrLoss), mmd = as.array(mmdLoss))
        }

        stage3Epochs <- as.integer(jointEpochs)
        for (epoch in seq_len(stage3Epochs)) {
            idx <- sample(nTrain)
            xShuffled <- xTrain[idx, , drop = FALSE]
            yShuffledOneHot <- yTrainOneHot[idx, , drop = FALSE]
            yShuffledLabels <- yTrain[idx]

            epochLosses <- list(total = 0, recon = 0, sw = 0, cluster = 0, contrastive = 0)
            for (batch in seq_len(nBatches)) {
                startIdx <- (batch - 1) * batchSize + 1
                endIdx <- min(batch * batchSize, nTrain)
                xBatch <- xShuffled[startIdx:endIdx, , drop = FALSE]
                yBatchOneHot <- yShuffledOneHot[startIdx:endIdx, , drop = FALSE]
                yBatchLabels <- yShuffledLabels[startIdx:endIdx]

                losses <- jointTrainStep(xBatch, yBatchOneHot, yBatchLabels)
                epochLosses$total <- epochLosses$total + losses$total
                epochLosses$recon <- epochLosses$recon + losses$recon
                epochLosses$sw <- epochLosses$sw + losses$sw
                epochLosses$cluster <- epochLosses$cluster + losses$cluster
                epochLosses$contrastive <- epochLosses$contrastive + losses$contrastive
            }

            if (epoch %% 10 == 0 || epoch == 1) {
                if (contrastiveLoss) {
                    message(sprintf("Stage 3 - Epoch %d/%d - total: %.4f, recon: %.4f, sw: %.4f, cluster: %.4f, contr: %.4f",
                                   epoch, stage3Epochs,
                                   epochLosses$total / nBatches,
                                   epochLosses$recon / nBatches,
                                   epochLosses$sw / nBatches,
                                   epochLosses$cluster / nBatches,
                                   epochLosses$contrastive / nBatches))
                } else {
                    message(sprintf("Stage 3 - Epoch %d/%d - total: %.4f, recon: %.4f, sw: %.4f, cluster: %.4f",
                                   epoch, stage3Epochs,
                                   epochLosses$total / nBatches,
                                   epochLosses$recon / nBatches,
                                   epochLosses$sw / nBatches,
                                   epochLosses$cluster / nBatches))
                }
            }
        }
    }

    message("\nTraining complete!")

    result <- list(
        encoder = encoder,
        decoder = decoder,
        clusterHead = clusterHead,
        history = history,
        nClusters = nClusters,
        latentDim = latentDim,
        supervised = supervised
    )

    if (useSwitchTab && !is.null(salientEncoder)) {
        result$salientEncoder <- salientEncoder
        result$mutualEncoder <- mutualEncoder
        result$switchDecoder <- switchDecoder
    }

    if (useAuxiliaryTask && !is.null(auxiliaryHead)) {
        result$auxiliaryHead <- auxiliaryHead
    }

    return(result)
}


#' Train SWAE-Clustering (Fast Version)
#'
#' A faster alternative to \code{\link{trainSWAEClustering}} that eliminates
#' per-batch R-to-TF tensor conversions. Uses the same R-level training loop
#' and shuffling as the original but pre-converts all data to TF tensors once
#' and uses TF gather for batching. Supports the core SWAE + DSW + supervised
#' clustering pipeline.
#'
#' @param trainData A data frame or matrix containing the training data.
#' @param useMarkers Character vector of column names to use as input features.
#' @param nClusters Integer, number of clusters (default 10).
#' @param pseudoLabels Optional vector of pseudo-labels for supervised training.
#' @param epochs Integer, total number of training epochs (default 100).
#' @param latentDim Integer, dimensionality of the latent space (default 16).
#' @param lambda Numeric, weight for Sliced Wasserstein distance loss (default 0.1).
#' @param lambdaCluster Numeric, weight for clustering loss (default 1.0).
#' @param numProjections Integer, number of random projections for SW (default 50).
#' @param batchSize Integer, training batch size (default 32).
#' @param learningRate Numeric, learning rate for optimizer (default 1e-3).
#' @param hiddenSizes Integer vector of length 2, hidden layer sizes (default c(23, 19)).
#' @param clusterWarmup Integer, epochs for Stage 1 SWAE warmup (default 10).
#' @param jointFineTune Logical, whether to do Stage 3 joint fine-tuning (default FALSE).
#' @param jointEpochs Integer, number of Stage 3 epochs (default 30).
#' @param labelSmoothing Numeric, label smoothing factor (default 0.0).
#' @param dropout Numeric, dropout rate (default 0.0).
#' @param useBatchNorm Logical, use batch normalization (default FALSE).
#' @param useLayerNorm Logical, use layer normalization (default FALSE).
#' @param weightDecay Numeric, weight decay for AdamW (default 0.0).
#' @param useDSW Logical, use Distributional Sliced Wasserstein (default FALSE).
#' @param dswMaxIter Integer, DSW inner optimization steps (default 10).
#' @param dswLambda Numeric, DSW diversity regularization strength (default 1.0).
#' @param verbose Integer, verbosity level (default 1).
#' @param seed Integer, random seed (default 1994).
#' @return A list containing trained encoder, decoder, clusterHead models, and training history.
#' @export
#' @examples
#' \dontrun{
#' model <- trainSWAEClusteringFast(trainData, useMarkers, nClusters = 11,
#'   pseudoLabels = labels, epochs = 100, latentDim = 16,
#'   hiddenSizes = c(23, 19), clusterWarmup = 20,
#'   useDSW = TRUE, dswMaxIter = 10, dswLambda = 1.0)
#' }
trainSWAEClusteringFast <- function(trainData, useMarkers, nClusters = 10L,
                                     pseudoLabels = NULL, epochs = 100L, latentDim = 16L,
                                     lambda = 0.1, lambdaCluster = 1.0, numProjections = 50L,
                                     batchSize = 32L, learningRate = 1e-3,
                                     hiddenSizes = c(23L, 19L), clusterWarmup = 10L,
                                     jointFineTune = FALSE, jointEpochs = 30L,
                                     labelSmoothing = 0.0,
                                     dropout = 0.0, useBatchNorm = FALSE, useLayerNorm = FALSE,
                                     weightDecay = 0.0,
                                     useDSW = FALSE, dswMaxIter = 10L, dswLambda = 1.0,
                                     verbose = 1L, seed = 1994) {

    tensorflow::set_random_seed(seed)
    tf <- tensorflow::tf
    set.seed(seed)

    # ========== DATA PREPARATION ==========
    xAll <- as.matrix(trainData[, useMarkers])
    inputDim <- ncol(xAll)
    nAll <- nrow(xAll)
    xTrain <- xAll
    nTrain <- nAll

    # Ensure integer types
    nClusters <- as.integer(nClusters)
    latentDim <- as.integer(latentDim)
    batchSize <- as.integer(batchSize)
    numProjections <- as.integer(numProjections)
    hiddenSizes <- as.integer(hiddenSizes)
    clusterWarmup <- as.integer(clusterWarmup)

    intermediateDim <- hiddenSizes[1]
    intermediateDim2 <- hiddenSizes[2]

    # Handle pseudo-labels
    supervised <- !is.null(pseudoLabels)
    yTrainOneHot <- NULL
    if (supervised) {
        if (!is.factor(pseudoLabels)) {
            pseudoLabels <- as.factor(pseudoLabels)
        }
        yTrain <- as.integer(pseudoLabels) - 1L
        if (any(is.na(yTrain))) stop("Invalid pseudo-labels: could not convert to integers")
        if (max(yTrain) >= nClusters) {
            stop(sprintf("Pseudo-labels have %d unique values but nClusters=%d",
                         max(yTrain) + 1, nClusters))
        }
        yTrainOneHot <- matrix(0, nrow = nTrain, ncol = nClusters)
        for (i in seq_len(nTrain)) {
            yTrainOneHot[i, yTrain[i] + 1] <- 1
        }
        message(sprintf("  Supervised mode: %d samples, %d clusters, label range %d-%d",
                         nTrain, nClusters, min(yTrain), max(yTrain)))
    }

    message("SWAE-Clustering Fast Architecture:")
    message(sprintf("  Input dim: %d", inputDim))
    message(sprintf("  Hidden: %d -> %d -> latent(%d)", intermediateDim, intermediateDim2, latentDim))
    message(sprintf("  Clusters: %d", nClusters))
    message(sprintf("  Training: %s mode", ifelse(supervised, "supervised", "self-supervised")))
    message(sprintf("  Lambda (SW): %.2f, Lambda (cluster): %.2f", lambda, lambdaCluster))
    message(sprintf("  Cluster warmup: %d epochs", clusterWarmup))
    if (useDSW) {
        message(sprintf("  Distributional SW: max_iter=%d, lambda=%.2f", dswMaxIter, dswLambda))
    }
    message("  Fast mode: compiled fit() with R shuffle + tf.data.Dataset (no TF shuffle)")

    # ========== BUILD ENCODER ==========
    encoderInput <- keras3::layer_input(shape = inputDim, name = "encoder_input")
    h <- encoderInput |>
        keras3::layer_dense(intermediateDim, activation = "relu", name = "enc_h1")
    if (useBatchNorm) h <- h |> keras3::layer_batch_normalization(name = "enc_bn1")
    if (useLayerNorm) h <- h |> keras3::layer_layer_normalization(name = "enc_ln1")
    if (dropout > 0) h <- h |> keras3::layer_dropout(rate = dropout, name = "enc_drop1")
    h <- h |>
        keras3::layer_dense(intermediateDim2, activation = "relu", name = "enc_h2")
    if (useBatchNorm) h <- h |> keras3::layer_batch_normalization(name = "enc_bn2")
    if (useLayerNorm) h <- h |> keras3::layer_layer_normalization(name = "enc_ln2")
    if (dropout > 0) h <- h |> keras3::layer_dropout(rate = dropout, name = "enc_drop2")
    h <- h |>
        keras3::layer_dense(latentDim, activation = "gelu", name = "latent")
    encoder <- keras3::keras_model(encoderInput, h, name = "encoder")

    # ========== BUILD DECODER ==========
    decoderInput <- keras3::layer_input(shape = latentDim, name = "decoder_input")
    decoded <- decoderInput |>
        keras3::layer_dense(intermediateDim2, activation = "relu", name = "dec_h1")
    if (useBatchNorm) decoded <- decoded |> keras3::layer_batch_normalization(name = "dec_bn1")
    if (useLayerNorm) decoded <- decoded |> keras3::layer_layer_normalization(name = "dec_ln1")
    if (dropout > 0) decoded <- decoded |> keras3::layer_dropout(rate = dropout, name = "dec_drop1")
    decoded <- decoded |>
        keras3::layer_dense(intermediateDim, activation = "relu", name = "dec_h2")
    if (useBatchNorm) decoded <- decoded |> keras3::layer_batch_normalization(name = "dec_bn2")
    if (useLayerNorm) decoded <- decoded |> keras3::layer_layer_normalization(name = "dec_ln2")
    if (dropout > 0) decoded <- decoded |> keras3::layer_dropout(rate = dropout, name = "dec_drop2")
    decoded <- decoded |>
        keras3::layer_dense(inputDim, activation = "sigmoid", name = "reconstruction")
    decoder <- keras3::keras_model(decoderInput, decoded, name = "decoder")

    # ========== BUILD CLUSTERING HEAD ==========
    clusterInput <- keras3::layer_input(shape = latentDim, name = "cluster_input")
    clusterLogits <- clusterInput |>
        keras3::layer_dense(32L, activation = "relu", name = "cluster_h1")
    if (dropout > 0) {
        clusterLogits <- clusterLogits |> keras3::layer_dropout(rate = dropout / 2, name = "cluster_drop")
    }
    clusterLogits <- clusterLogits |>
        keras3::layer_dense(nClusters, activation = "softmax", name = "cluster_output")
    clusterHead <- keras3::keras_model(clusterInput, clusterLogits, name = "cluster_head")

    # ========== BUILD DSW PROJECTION NETWORK (optional) ==========
    dswProjectionNet <- NULL
    dswProjectionOptimizer <- NULL
    if (useDSW) {
        dswInput <- keras3::layer_input(shape = latentDim, name = "dsw_input")
        dswH <- dswInput |>
            keras3::layer_dense(latentDim, activation = "relu", name = "dsw_h1") |>
            keras3::layer_dense(latentDim, activation = "linear", name = "dsw_output")
        dswProjectionNet <- keras3::keras_model(dswInput, dswH, name = "dsw_projection")
        dswProjectionOptimizer <- keras3::optimizer_adam(learning_rate = 1e-4)
        message("  DSW projection network built")
    }

    # ========== SETUP OPTIMIZERS ==========
    if (weightDecay > 0) {
        optimizerAE <- keras3::optimizer_adam_w(learning_rate = learningRate, weight_decay = weightDecay)
    } else {
        optimizerAE <- keras3::optimizer_adam(learning_rate = learningRate)
    }

    # Pre-compute TF constant for inputDim (used in compiled train_step)
    inputDimFloat <- tf$constant(as.numeric(inputDim), dtype = tf$float32)

    # ========== DSW INNER LOOP IN PURE PYTHON (avoids R bridge overhead) ==========
    dswInnerOptimizePy <- NULL
    if (useDSW) {
        # Set projection net/optimizer as Python globals for closure capture
        # NOTE: must use local py variable — reticulate::py$name <- val fails
        #       because R's $<- can't write back through :: accessor
        py <- reticulate::py
        reticulate::py_run_string("import tensorflow as tf")
        py$`_dsw_proj_net` <- dswProjectionNet
        py$`_dsw_proj_opt` <- dswProjectionOptimizer
        py$`_dsw_max_iter` <- as.integer(dswMaxIter)
        py$`_dsw_lam` <- dswLambda
        py$`_dsw_num_proj` <- as.integer(numProjections)

        reticulate::py_run_string("
def _dsw_inner_optimize(pro, z_detach, prior_detach):
    '''Run DSW projection optimization eagerly in Python (no R bridge overhead).
    Called via tf.py_function from compiled train_step.'''
    for i in range(_dsw_max_iter):
        with tf.GradientTape() as tape:
            projections = _dsw_proj_net(pro, training=True)
            projections = projections / tf.norm(projections, axis=1, keepdims=True)
            cosine = tf.matmul(projections, projections, transpose_b=True)
            identity = tf.eye(_dsw_num_proj)
            reg = _dsw_lam * tf.reduce_mean(tf.square(cosine - identity))
            z_proj = tf.matmul(z_detach, projections, transpose_b=True)
            prior_proj = tf.matmul(prior_detach, projections, transpose_b=True)
            z_sorted = tf.sort(z_proj, axis=0)
            prior_sorted = tf.sort(prior_proj, axis=0)
            wd = tf.reduce_mean(tf.square(z_sorted - prior_sorted))
            inner_loss = reg - wd
        grads = tape.gradient(inner_loss, _dsw_proj_net.trainable_variables)
        _dsw_proj_opt.apply_gradients(zip(grads, _dsw_proj_net.trainable_variables))
    return tf.constant(0.0)
")
        dswInnerOptimizePy <- py$`_dsw_inner_optimize`
        message("  DSW inner loop: Python py_function (eager within compiled graph)")
    }

    # ========== STAGE 1: SWAE MODEL CLASS (compiled, DSW via py_function) ==========
    SWAEModel <- keras3::new_model_class(
        classname = "SWAEModel",
        initialize = function(enc, dec, ...) {
            super$initialize(...)
            self$enc <- enc
            self$dec <- dec
            self$total_loss_tracker <- keras3::metric_mean(name = "total_loss")
            self$recon_loss_tracker <- keras3::metric_mean(name = "recon_loss")
            self$sw_loss_tracker <- keras3::metric_mean(name = "sw_loss")
        },
        metrics = keras3::mark_active(function() {
            list(self$total_loss_tracker, self$recon_loss_tracker, self$sw_loss_tracker)
        }),
        train_step = function(data) {
            with(tf$GradientTape() %as% tape, {
                z <- self$enc(data, training = TRUE)
                xRecon <- self$dec(z, training = TRUE)
                reconLoss <- tf$reduce_mean(tf$square(data - xRecon)) * inputDimFloat

                priorSamples <- tf$random$normal(shape = tf$shape(z))
                if (useDSW) {
                    # Generate random projections (use known integer dims to avoid shape issues)
                    pro <- tf$random$normal(shape = c(as.integer(numProjections), as.integer(latentDim)))
                    pro <- pro / tf$norm(pro, axis = 1L, keepdims = TRUE)

                    # Run DSW inner optimization eagerly via py_function
                    zDetach <- tf$stop_gradient(z)
                    priorDetach <- tf$stop_gradient(priorSamples)
                    dummy <- tf$py_function(
                        dswInnerOptimizePy,
                        list(pro, zDetach, priorDetach),
                        tf$float32)
                    # py_function outputs have unknown shape; set it to scalar
                    dummy <- tf$ensure_shape(dummy, shape = list())

                    # Force dependency: ensure inner opt completes before reading weights
                    proDep <- pro + dummy * tf$constant(0.0)

                    # Final DSW computation in graph mode (differentiable w.r.t. z)
                    projections <- dswProjectionNet(proDep, training = FALSE)
                    projections <- projections / tf$norm(projections, axis = 1L, keepdims = TRUE)
                    zProj <- tf$matmul(z, projections, transpose_b = TRUE)
                    priorProj <- tf$matmul(priorSamples, projections, transpose_b = TRUE)
                    zSorted <- tf$sort(zProj, axis = 0L)
                    priorSorted <- tf$sort(priorProj, axis = 0L)
                    swLoss <- tf$reduce_mean(tf$square(zSorted - priorSorted))
                } else {
                    swLoss <- slicedWassersteinDistance(z, priorSamples, numProjections)
                }

                totalLoss <- reconLoss + lambda * swLoss
            })

            grads <- tape$gradient(totalLoss, self$trainable_weights)
            self$optimizer$apply_gradients(keras3::zip_lists(grads, self$trainable_weights))

            self$total_loss_tracker$update_state(totalLoss)
            self$recon_loss_tracker$update_state(reconLoss)
            self$sw_loss_tracker$update_state(swLoss)
            list(total_loss = self$total_loss_tracker$result(),
                 recon_loss = self$recon_loss_tracker$result(),
                 sw_loss = self$sw_loss_tracker$result())
        }
    )

    swaeModel <- SWAEModel(encoder, decoder)
    keras3::compile(swaeModel, optimizer = optimizerAE)

    # ========== STAGE 1: TRAIN SWAE (compiled fit, R shuffle, DSW via py_function) ==========
    message(sprintf("\n=== STAGE 1: Training SWAE (%d epochs, compiled+py_function) ===", clusterWarmup))

    nBatches <- ceiling(nTrain / batchSize)
    for (epoch in seq_len(clusterWarmup)) {
        idx <- sample(nTrain)
        xShuffled <- xTrain[idx, , drop = FALSE]
        xTensor <- tf$constant(xShuffled, dtype = tf$float32)
        dataset <- tf$data$Dataset$from_tensor_slices(xTensor)$batch(batchSize)$prefetch(tf$data$AUTOTUNE)
        h <- swaeModel$fit(dataset, epochs = 1L, verbose = 0L)

        if (epoch %% 10 == 0 || epoch == 1) {
            message(sprintf("Stage 1 - Epoch %d/%d - loss: %.4f",
                            epoch, clusterWarmup, as.numeric(h$history$total_loss)))
        }
    }

    # ========== STAGE 2: CLUSTER HEAD MODEL CLASS (compiled) ==========
    # Ensure yTrainOneHot exists (dummy for unsupervised DEC loss)
    if (!supervised) {
        yTrainOneHot <- matrix(0, nrow = nTrain, ncol = nClusters)
    }

    ClusterModel <- keras3::new_model_class(
        classname = "ClusterModel",
        initialize = function(head, ...) {
            super$initialize(...)
            self$head <- head
            self$cluster_loss_tracker <- keras3::metric_mean(name = "cluster_loss")
        },
        metrics = keras3::mark_active(function() {
            list(self$cluster_loss_tracker)
        }),
        train_step = function(data) {
            x <- data[[1]]
            yOneHot <- data[[2]]
            z <- tf$stop_gradient(encoder(x, training = FALSE))

            with(tf$GradientTape() %as% tape, {
                probs <- self$head(z, training = TRUE)
                if (supervised) {
                    smoothed <- yOneHot * (1.0 - labelSmoothing) +
                        labelSmoothing / tf$cast(nClusters, tf$float32)
                    clusterLoss <- tf$reduce_mean(
                        tf$keras$losses$categorical_crossentropy(smoothed, probs))
                } else {
                    q <- probs + 1e-10
                    f <- tf$reduce_sum(q, axis = 0L, keepdims = TRUE)
                    p <- tf$square(q) / f
                    p <- p / tf$reduce_sum(p, axis = 1L, keepdims = TRUE)
                    clusterLoss <- tf$reduce_mean(
                        tf$reduce_sum(p * tf$math$log(p / q), axis = 1L))
                }
            })

            grads <- tape$gradient(clusterLoss, self$trainable_weights)
            self$optimizer$apply_gradients(keras3::zip_lists(grads, self$trainable_weights))
            self$cluster_loss_tracker$update_state(clusterLoss)
            list(cluster_loss = self$cluster_loss_tracker$result())
        }
    )

    optimizerCluster <- keras3::optimizer_adam(learning_rate = learningRate)
    clusterModel <- ClusterModel(clusterHead)
    keras3::compile(clusterModel, optimizer = optimizerCluster)

    # ========== STAGE 2: TRAIN CLUSTERING HEAD (R shuffle + compiled fit) ==========
    stage2Epochs <- as.integer(epochs - clusterWarmup)
    message(sprintf("\n=== STAGE 2: Training clustering head (%d epochs, compiled) ===", stage2Epochs))

    history <- list(cluster = c())
    for (epoch in seq_len(stage2Epochs)) {
        idx <- sample(nTrain)
        xShuffled <- xTrain[idx, , drop = FALSE]
        yShuffled <- yTrainOneHot[idx, , drop = FALSE]
        xTensor <- tf$constant(xShuffled, dtype = tf$float32)
        yTensor <- tf$constant(yShuffled, dtype = tf$float32)

        dataset <- tf$data$Dataset$from_tensor_slices(
            reticulate::tuple(xTensor, yTensor))$batch(batchSize)$prefetch(tf$data$AUTOTUNE)
        h <- clusterModel$fit(dataset, epochs = 1L, verbose = 0L)

        clusterLossVal <- as.numeric(h$history$cluster_loss)
        history$cluster <- c(history$cluster, clusterLossVal)

        if (epoch %% 10 == 0 || epoch == 1) {
            message(sprintf("Stage 2 - Epoch %d/%d - cluster_loss: %.4f",
                            epoch, stage2Epochs, clusterLossVal))
        }
    }

    # ========== STAGE 3: JOINT FINE-TUNING (optional, compiled) ==========
    if (supervised && jointFineTune) {
        JointModel <- keras3::new_model_class(
            classname = "JointModel",
            initialize = function(enc, dec, head, ...) {
                super$initialize(...)
                self$enc <- enc
                self$dec <- dec
                self$head <- head
                self$total_loss_tracker <- keras3::metric_mean(name = "total_loss")
            },
            metrics = keras3::mark_active(function() {
                list(self$total_loss_tracker)
            }),
            train_step = function(data) {
                x <- data[[1]]
                yOneHot <- data[[2]]

                with(tf$GradientTape() %as% tape, {
                    z <- self$enc(x, training = TRUE)
                    xRecon <- self$dec(z, training = TRUE)
                    probs <- self$head(z, training = TRUE)

                    reconLoss <- tf$reduce_mean(tf$square(x - xRecon)) * inputDimFloat

                    priorSamples <- tf$random$normal(shape = tf$shape(z))
                    if (useDSW) {
                        # Same py_function approach as Stage 1
                        pro <- tf$random$normal(shape = c(as.integer(numProjections), as.integer(latentDim)))
                        pro <- pro / tf$norm(pro, axis = 1L, keepdims = TRUE)
                        zDetach <- tf$stop_gradient(z)
                        priorDetach <- tf$stop_gradient(priorSamples)
                        dummy <- tf$py_function(
                            dswInnerOptimizePy,
                            list(pro, zDetach, priorDetach),
                            tf$float32)
                        dummy <- tf$ensure_shape(dummy, shape = list())
                        proDep <- pro + dummy * tf$constant(0.0)
                        projections <- dswProjectionNet(proDep, training = FALSE)
                        projections <- projections / tf$norm(projections, axis = 1L, keepdims = TRUE)
                        zProj <- tf$matmul(z, projections, transpose_b = TRUE)
                        priorProj <- tf$matmul(priorSamples, projections, transpose_b = TRUE)
                        zSorted <- tf$sort(zProj, axis = 0L)
                        priorSorted <- tf$sort(priorProj, axis = 0L)
                        swLoss <- tf$reduce_mean(tf$square(zSorted - priorSorted))
                    } else {
                        swLoss <- slicedWassersteinDistance(z, priorSamples, numProjections)
                    }

                    smoothedY <- yOneHot * (1.0 - labelSmoothing) +
                        labelSmoothing / tf$cast(as.integer(nClusters), tf$float32)
                    clusterLoss <- tf$reduce_mean(
                        tf$keras$losses$categorical_crossentropy(smoothedY, probs))
                    totalLoss <- reconLoss + lambda * swLoss + lambdaCluster * clusterLoss
                })

                grads <- tape$gradient(totalLoss, self$trainable_weights)
                self$optimizer$apply_gradients(keras3::zip_lists(grads, self$trainable_weights))
                self$total_loss_tracker$update_state(totalLoss)
                list(total_loss = self$total_loss_tracker$result())
            }
        )

        message(sprintf("\n=== STAGE 3: Joint fine-tuning (%d epochs, compiled) ===", jointEpochs))
        optimizerJoint <- keras3::optimizer_adam(learning_rate = learningRate * 0.1, clipnorm = 1.0)
        jointModel <- JointModel(encoder, decoder, clusterHead)
        keras3::compile(jointModel, optimizer = optimizerJoint)

        stage3Epochs <- as.integer(jointEpochs)
        for (epoch in seq_len(stage3Epochs)) {
            idx <- sample(nTrain)
            xShuffled <- xTrain[idx, , drop = FALSE]
            yShuffled <- yTrainOneHot[idx, , drop = FALSE]
            xTensor <- tf$constant(xShuffled, dtype = tf$float32)
            yTensor <- tf$constant(yShuffled, dtype = tf$float32)

            dataset <- tf$data$Dataset$from_tensor_slices(
                reticulate::tuple(xTensor, yTensor))$batch(batchSize)$prefetch(tf$data$AUTOTUNE)
            h <- jointModel$fit(dataset, epochs = 1L, verbose = 0L)

            if (epoch %% 10 == 0 || epoch == 1) {
                message(sprintf("Stage 3 - Epoch %d/%d - total_loss: %.4f",
                                epoch, stage3Epochs, as.numeric(h$history$total_loss)))
            }
        }
    }

    message("\nTraining complete!")

    result <- list(
        encoder = encoder,
        decoder = decoder,
        clusterHead = clusterHead,
        history = history,
        nClusters = nClusters,
        latentDim = latentDim,
        supervised = supervised
    )

    return(result)
}


#' Train VQ-VAE with Clustering
#'
#' Trains a Vector Quantized Variational Autoencoder where codebook entries
#' serve as learned cell type prototypes. Compatible with encodeSWAE/decodeSWAE.
#'
#' @param trainData Data frame containing training data.
#' @param useMarkers Character vector of marker column names.
#' @param nClusters Integer, number of codebook entries / clusters.
#' @param pseudoLabels Optional factor of pseudo-labels for supervised fine-tuning.
#' @param epochs Integer, total training epochs.
#' @param latentDim Integer, latent / codebook vector dimension.
#' @param beta Numeric, commitment loss weight (default 0.25).
#' @param lambdaCluster Numeric, supervised cluster loss weight.
#' @param batchSize Integer, training batch size.
#' @param learningRate Numeric, Adam learning rate.
#' @param hiddenSizes Integer vector of length 2 for hidden layer sizes.
#' @param clusterWarmup Integer, unsupervised VQ-VAE epochs before adding cluster loss.
#' @param labelSmoothing Numeric, label smoothing for supervised loss.
#' @param dropout Numeric, dropout rate.
#' @param useBatchNorm Logical, use batch normalization.
#' @param useLayerNorm Logical, use layer normalization.
#' @param weightDecay Numeric, weight decay for AdamW.
#' @param codebookReset Logical, reset unused codebook entries periodically.
#' @param resetInterval Integer, how often (epochs) to check for dead codes.
#' @param verbose Integer, verbosity level.
#' @param seed Integer, random seed.
#' @return List with encoder, decoder, clusterHead, codebook, history, metadata.
#' @export
trainVQVAEClustering <- function(trainData, useMarkers, nClusters = 10L,
                                  pseudoLabels = NULL, epochs = 100L, latentDim = 16L,
                                  beta = 0.25, lambdaCluster = 1.0,
                                  batchSize = 64L, learningRate = 1e-3,
                                  hiddenSizes = c(23L, 19L), clusterWarmup = 20L,
                                  labelSmoothing = 0.0,
                                  dropout = 0.0, useBatchNorm = FALSE, useLayerNorm = FALSE,
                                  weightDecay = 0.0,
                                  codebookReset = TRUE, resetInterval = 10L,
                                  softQuantize = FALSE, temperature = 1.0,
                                  skipVQDecode = FALSE,
                                  useEMA = FALSE, emaDecay = 0.99,
                                  useSkipConnections = FALSE,
                                  reconWeights = NULL,
                                  adversarial = FALSE, adversarialLambda = 1.0,
                                  sampleLabels = NULL,
                                  verbose = 1L, seed = 1994) {

    tensorflow::set_random_seed(seed)
    tf <- tensorflow::tf
    set.seed(seed)

    # ========== DATA PREPARATION ==========
    xAll <- as.matrix(trainData[, useMarkers])
    inputDim <- ncol(xAll)
    nAll <- nrow(xAll)
    xTrain <- xAll
    nTrain <- nAll

    nClusters <- as.integer(nClusters)
    latentDim <- as.integer(latentDim)
    batchSize <- as.integer(batchSize)
    hiddenSizes <- as.integer(hiddenSizes)
    clusterWarmup <- as.integer(clusterWarmup)

    intermediateDim <- hiddenSizes[1]
    intermediateDim2 <- hiddenSizes[2]

    supervised <- !is.null(pseudoLabels)
    yTrainOneHot <- NULL
    if (supervised) {
        if (!is.factor(pseudoLabels)) pseudoLabels <- as.factor(pseudoLabels)
        yTrain <- as.integer(pseudoLabels) - 1L
        if (any(is.na(yTrain))) stop("Invalid pseudo-labels")
        if (max(yTrain) >= nClusters) {
            stop(sprintf("Pseudo-labels have %d unique values but nClusters=%d",
                         max(yTrain) + 1, nClusters))
        }
        yTrainOneHot <- matrix(0, nrow = nTrain, ncol = nClusters)
        for (i in seq_len(nTrain)) yTrainOneHot[i, yTrain[i] + 1] <- 1
        message(sprintf("  Supervised mode: %d samples, %d clusters", nTrain, nClusters))
    }

    # ========== WEIGHTED RECONSTRUCTION ==========
    reconWeightsTF <- NULL
    if (!is.null(reconWeights)) {
        reconWeights <- as.numeric(reconWeights)
        if (length(reconWeights) != inputDim) stop("reconWeights must have length == number of markers")
        reconWeights <- reconWeights / mean(reconWeights)  # normalize so mean weight = 1
        reconWeightsTF <- tf$constant(matrix(reconWeights, nrow = 1), dtype = tf$float32)
    }

    # ========== ADVERSARIAL BATCH REMOVAL ==========
    nSamplesAdv <- 0L
    sampleLabelsTrain <- NULL
    if (adversarial) {
        if (is.null(sampleLabels)) stop("sampleLabels required when adversarial=TRUE")
        sampleLabelsFactor <- as.factor(sampleLabels)
        nSamplesAdv <- as.integer(nlevels(sampleLabelsFactor))
        sampleLabelsTrain <- as.integer(sampleLabelsFactor) - 1L
        advLambdaFloat <- tf$constant(adversarialLambda, dtype = tf$float32)
        nSamplesAdvFloat <- tf$constant(as.numeric(nSamplesAdv), dtype = tf$float32)
    }

    message("VQ-VAE Clustering Architecture:")
    message(sprintf("  Input dim: %d", inputDim))
    message(sprintf("  Hidden: %d -> %d -> latent(%d)", intermediateDim, intermediateDim2, latentDim))
    message(sprintf("  Codebook: %d entries x %d dims", nClusters, latentDim))
    message(sprintf("  Beta (commitment): %.2f", beta))
    message(sprintf("  Soft quantize: %s (temp=%.2f), Skip VQ decode: %s, EMA: %s (decay=%.3f)",
                    softQuantize, temperature, skipVQDecode, useEMA, emaDecay))
    if (useSkipConnections) message("  Skip connections: enabled")
    if (!is.null(reconWeights)) message("  Weighted reconstruction: enabled")
    if (adversarial) message(sprintf("  Adversarial: %d samples, lambda=%.2f", nSamplesAdv, adversarialLambda))
    if (supervised) {
        message(sprintf("  Supervised: warmup=%d epochs, then +cluster loss (lambda=%.2f)",
                        clusterWarmup, lambdaCluster))
    }

    # ========== BUILD ENCODER ==========
    encoderInput <- keras3::layer_input(shape = inputDim, name = "encoder_input")
    h <- encoderInput |>
        keras3::layer_dense(intermediateDim, activation = "relu", name = "enc_h1")
    if (useBatchNorm) h <- h |> keras3::layer_batch_normalization(name = "enc_bn1")
    if (useLayerNorm) h <- h |> keras3::layer_layer_normalization(name = "enc_ln1")
    if (dropout > 0) h <- h |> keras3::layer_dropout(rate = dropout, name = "enc_drop1")
    if (useSkipConnections) {
        skip1 <- keras3::layer_dense(encoderInput, units = intermediateDim,
                                      activation = NULL, use_bias = FALSE, name = "enc_skip1")
        h <- keras3::layer_add(list(h, skip1))
    }
    h <- h |>
        keras3::layer_dense(intermediateDim2, activation = "relu", name = "enc_h2")
    if (useBatchNorm) h <- h |> keras3::layer_batch_normalization(name = "enc_bn2")
    if (useLayerNorm) h <- h |> keras3::layer_layer_normalization(name = "enc_ln2")
    if (dropout > 0) h <- h |> keras3::layer_dropout(rate = dropout, name = "enc_drop2")
    if (useSkipConnections) {
        skip2_in <- encoderInput |>
            keras3::layer_dense(intermediateDim2, activation = NULL, use_bias = FALSE, name = "enc_skip2")
        h <- keras3::layer_add(list(h, skip2_in))
    }
    h <- h |>
        keras3::layer_dense(latentDim, activation = "gelu", name = "latent")
    encoder <- keras3::keras_model(encoderInput, h, name = "encoder")

    # ========== BUILD DECODER ==========
    decoderInput <- keras3::layer_input(shape = latentDim, name = "decoder_input")
    decoded <- decoderInput |>
        keras3::layer_dense(intermediateDim2, activation = "relu", name = "dec_h1")
    if (useBatchNorm) decoded <- decoded |> keras3::layer_batch_normalization(name = "dec_bn1")
    if (useLayerNorm) decoded <- decoded |> keras3::layer_layer_normalization(name = "dec_ln1")
    if (dropout > 0) decoded <- decoded |> keras3::layer_dropout(rate = dropout, name = "dec_drop1")
    if (useSkipConnections) {
        dec_skip1 <- keras3::layer_dense(decoderInput, units = intermediateDim2,
                                          activation = NULL, use_bias = FALSE, name = "dec_skip1")
        decoded <- keras3::layer_add(list(decoded, dec_skip1))
    }
    decoded <- decoded |>
        keras3::layer_dense(intermediateDim, activation = "relu", name = "dec_h2")
    if (useBatchNorm) decoded <- decoded |> keras3::layer_batch_normalization(name = "dec_bn2")
    if (useLayerNorm) decoded <- decoded |> keras3::layer_layer_normalization(name = "dec_ln2")
    if (dropout > 0) decoded <- decoded |> keras3::layer_dropout(rate = dropout, name = "dec_drop2")
    if (useSkipConnections) {
        dec_skip2 <- keras3::layer_dense(decoderInput, units = intermediateDim,
                                          activation = NULL, use_bias = FALSE, name = "dec_skip2")
        decoded <- keras3::layer_add(list(decoded, dec_skip2))
    }
    decoded <- decoded |>
        keras3::layer_dense(inputDim, activation = "sigmoid", name = "reconstruction")
    decoder <- keras3::keras_model(decoderInput, decoded, name = "decoder")

    # ========== BUILD ADVERSARIAL DISCRIMINATOR ==========
    discriminator <- NULL
    optimizerDisc <- NULL
    if (adversarial) {
        discInput <- keras3::layer_input(shape = latentDim, name = "disc_input")
        dh <- discInput |>
            keras3::layer_dense(32L, activation = "relu", name = "disc_h1") |>
            keras3::layer_dense(32L, activation = "relu", name = "disc_h2") |>
            keras3::layer_dense(nSamplesAdv, activation = "softmax", name = "disc_out")
        discriminator <- keras3::keras_model(discInput, dh, name = "discriminator")
        optimizerDisc <- keras3::optimizer_adam(learning_rate = learningRate)
    }

    # ========== SETUP OPTIMIZER ==========
    if (weightDecay > 0) {
        optimizerAE <- keras3::optimizer_adam_w(learning_rate = learningRate, weight_decay = weightDecay)
    } else {
        optimizerAE <- keras3::optimizer_adam(learning_rate = learningRate)
    }

    inputDimFloat <- tf$constant(as.numeric(inputDim), dtype = tf$float32)
    betaFloat <- tf$constant(beta, dtype = tf$float32)
    temperatureFloat <- tf$constant(temperature, dtype = tf$float32)
    emaDecayFloat <- tf$constant(emaDecay, dtype = tf$float32)
    nKInt <- as.integer(nClusters)

    # ========== VQ-VAE MODEL CLASS (Stage 1: unsupervised) ==========
    VQVAEModel <- keras3::new_model_class(
        classname = "VQVAEModel",
        initialize = function(enc, dec, nK, latDim, ...) {
            super$initialize(...)
            self$enc <- enc
            self$dec <- dec
            # Codebook: trainable only when NOT using EMA
            self$vq_codebook <- self$add_weight(
                name = "vq_codebook",
                shape = as.integer(c(nK, latDim)),
                initializer = keras3::initializer_random_uniform(
                    minval = -1.0 / nK, maxval = 1.0 / nK),
                trainable = !useEMA
            )
            if (useEMA) {
                self$ema_count <- self$add_weight(
                    name = "ema_count", shape = list(as.integer(nK)),
                    initializer = keras3::initializer_ones(), trainable = FALSE)
                self$ema_sum <- self$add_weight(
                    name = "ema_sum", shape = as.integer(c(nK, latDim)),
                    initializer = keras3::initializer_zeros(), trainable = FALSE)
            }
            self$total_loss_tracker <- keras3::metric_mean(name = "total_loss")
            self$recon_loss_tracker <- keras3::metric_mean(name = "recon_loss")
            self$commit_loss_tracker <- keras3::metric_mean(name = "commit_loss")
        },
        metrics = keras3::mark_active(function() {
            list(self$total_loss_tracker, self$recon_loss_tracker,
                 self$commit_loss_tracker)
        }),
        train_step = function(data) {
            with(tf$GradientTape() %as% tape, {
                z_e <- self$enc(data, training = TRUE)

                # VQ: L2 distances to codebook entries
                z_e_exp <- tf$expand_dims(z_e, axis = 1L)
                cb_exp <- tf$expand_dims(self$vq_codebook, axis = 0L)
                distances <- tf$reduce_sum(tf$square(z_e_exp - cb_exp), axis = 2L)

                if (softQuantize) {
                    # Soft quantization: weighted avg of codebook entries
                    weights <- tf$nn$softmax(-distances / temperatureFloat, axis = 1L)
                    z_q <- tf$linalg$matmul(weights, self$vq_codebook)
                } else {
                    # Hard quantization: nearest neighbor
                    indices <- tf$argmin(distances, axis = 1L)
                    z_q <- tf$gather(self$vq_codebook, indices)
                }

                if (skipVQDecode) {
                    # Decode from continuous z_e (no quantization loss)
                    x_recon <- self$dec(z_e, training = TRUE)
                } else {
                    # Straight-through estimator
                    z_q_st <- z_e + tf$stop_gradient(z_q - z_e)
                    x_recon <- self$dec(z_q_st, training = TRUE)
                }

                # Losses
                if (!is.null(reconWeightsTF)) {
                    reconLoss <- tf$reduce_mean(tf$square(data - x_recon) * reconWeightsTF) * inputDimFloat
                } else {
                    reconLoss <- tf$reduce_mean(tf$square(data - x_recon)) * inputDimFloat
                }
                commitLoss <- tf$reduce_mean(tf$square(z_e - tf$stop_gradient(z_q)))

                if (useEMA) {
                    totalLoss <- reconLoss + betaFloat * commitLoss
                } else {
                    codebookLoss <- tf$reduce_mean(tf$square(tf$stop_gradient(z_e) - z_q))
                    totalLoss <- reconLoss + betaFloat * commitLoss + codebookLoss
                }
            })

            grads <- tape$gradient(totalLoss, self$trainable_weights)
            self$optimizer$apply_gradients(keras3::zip_lists(grads, self$trainable_weights))

            # EMA codebook update (outside gradient tape)
            if (useEMA) {
                hardIdx <- tf$argmin(distances, axis = 1L)
                oneHotAssign <- tf$one_hot(hardIdx, depth = nKInt)
                n_k <- tf$reduce_sum(oneHotAssign, axis = 0L)
                sum_k <- tf$linalg$matmul(tf$transpose(oneHotAssign), z_e)
                self$ema_count$assign(emaDecayFloat * self$ema_count +
                                      (1.0 - emaDecayFloat) * n_k)
                self$ema_sum$assign(emaDecayFloat * self$ema_sum +
                                    (1.0 - emaDecayFloat) * sum_k)
                updated_cb <- self$ema_sum / tf$expand_dims(self$ema_count + 1e-5, axis = 1L)
                self$vq_codebook$assign(updated_cb)
            }

            self$total_loss_tracker$update_state(totalLoss)
            self$recon_loss_tracker$update_state(reconLoss)
            self$commit_loss_tracker$update_state(commitLoss)

            list(total_loss = self$total_loss_tracker$result(),
                 recon_loss = self$recon_loss_tracker$result(),
                 commit_loss = self$commit_loss_tracker$result())
        }
    )

    vqvaeModel <- VQVAEModel(enc = encoder, dec = decoder,
                              nK = nClusters, latDim = latentDim)
    keras3::compile(vqvaeModel, optimizer = optimizerAE)

    # Sync EMA state to initial codebook
    if (useEMA) {
        init_cb <- vqvaeModel$vq_codebook$numpy()
        vqvaeModel$ema_sum$assign(tf$constant(init_cb, dtype = tf$float32))
    }
    message(sprintf("  Codebook initialized: %d x %d (EMA=%s)", nClusters, latentDim, useEMA))

    # ========== STAGE 1: VQ-VAE TRAINING ==========
    if (supervised) {
        stage1Epochs <- clusterWarmup
    } else {
        stage1Epochs <- as.integer(epochs)
    }
    message(sprintf("\n=== STAGE 1: VQ-VAE training (%d epochs) ===", stage1Epochs))

    history <- list(total = c(), recon = c(), commit = c())
    for (epoch in seq_len(stage1Epochs)) {
        idx <- sample(nTrain)
        xShuffled <- xTrain[idx, , drop = FALSE]
        xTensor <- tf$constant(xShuffled, dtype = tf$float32)
        dataset <- tf$data$Dataset$from_tensor_slices(xTensor)$batch(batchSize)$prefetch(tf$data$AUTOTUNE)
        h <- vqvaeModel$fit(dataset, epochs = 1L, verbose = 0L)

        history$total <- c(history$total, as.numeric(h$history$total_loss))
        history$recon <- c(history$recon, as.numeric(h$history$recon_loss))

        # Dead code reset
        if (codebookReset && (epoch %% resetInterval == 0)) {
            z_all <- as.matrix(encoder(tf$constant(xTrain, dtype = tf$float32), training = FALSE))
            z_all_t <- tf$constant(z_all, dtype = tf$float32)
            z_exp <- tf$expand_dims(z_all_t, axis = 1L)
            cb_exp <- tf$expand_dims(vqvaeModel$vq_codebook, axis = 0L)
            dists <- tf$reduce_sum(tf$square(z_exp - cb_exp), axis = 2L)
            assignments <- as.integer(as.matrix(tf$argmin(dists, axis = 1L)))
            counts <- tabulate(assignments + 1L, nbins = nClusters)
            dead <- which(counts == 0)
            if (length(dead) > 0) {
                message(sprintf("  Epoch %d: resetting %d dead codebook entries", epoch, length(dead)))
                randIdx <- sample(nTrain, length(dead), replace = TRUE)
                newVals <- z_all[randIdx, , drop = FALSE] +
                    matrix(rnorm(length(dead) * latentDim, sd = 0.01),
                           nrow = length(dead))
                cbNumpy <- as.matrix(vqvaeModel$vq_codebook$numpy())
                cbNumpy[dead, ] <- newVals
                vqvaeModel$vq_codebook$assign(tf$constant(cbNumpy, dtype = tf$float32))
                if (useEMA) {
                    emaSumNp <- as.matrix(vqvaeModel$ema_sum$numpy())
                    emaCntNp <- as.numeric(vqvaeModel$ema_count$numpy())
                    emaSumNp[dead, ] <- newVals
                    emaCntNp[dead] <- 1.0
                    vqvaeModel$ema_sum$assign(tf$constant(emaSumNp, dtype = tf$float32))
                    vqvaeModel$ema_count$assign(tf$constant(emaCntNp, dtype = tf$float32))
                }
            }
        }

        # Adversarial batch removal step
        if (adversarial) {
            z_all_adv <- encoder(tf$constant(xTrain[idx, , drop = FALSE], dtype = tf$float32), training = FALSE)
            sLabels <- tf$constant(sampleLabelsTrain[idx], dtype = tf$int32)
            sOneHot <- tf$one_hot(sLabels, depth = as.integer(nSamplesAdv))

            # Step 1: Train discriminator to classify sample_id from z_e
            for (dStep in 1:3) {
                with(tf$GradientTape() %as% dtape, {
                    discPreds <- discriminator(tf$stop_gradient(z_all_adv), training = TRUE)
                    discLoss <- tf$reduce_mean(
                        tf$keras$losses$categorical_crossentropy(sOneHot, discPreds))
                })
                dgrads <- dtape$gradient(discLoss, discriminator$trainable_weights)
                optimizerDisc$apply_gradients(
                    keras3::zip_lists(dgrads, discriminator$trainable_weights))
            }

            # Step 2: Update encoder to FOOL discriminator (gradient reversal)
            with(tf$GradientTape() %as% etape, {
                z_adv <- encoder(tf$constant(xTrain[idx, , drop = FALSE], dtype = tf$float32), training = TRUE)
                discPreds2 <- discriminator(z_adv, training = FALSE)
                # Maximize discriminator loss = minimize negative entropy
                # Uniform prediction = maximum confusion
                uniformTarget <- tf$ones_like(discPreds2) / nSamplesAdvFloat
                advLoss <- -advLambdaFloat * tf$reduce_mean(
                    tf$keras$losses$categorical_crossentropy(uniformTarget, discPreds2))
            })
            egrads <- etape$gradient(advLoss, encoder$trainable_weights)
            # Filter out None gradients
            validGrads <- list()
            validWeights <- list()
            for (gi in seq_along(egrads)) {
                if (!is.null(egrads[[gi]])) {
                    validGrads[[length(validGrads) + 1]] <- egrads[[gi]]
                    validWeights[[length(validWeights) + 1]] <- encoder$trainable_weights[[gi]]
                }
            }
            if (length(validGrads) > 0) {
                optimizerAE$apply_gradients(keras3::zip_lists(validGrads, validWeights))
            }
        }

        if (epoch %% 10 == 0 || epoch == 1) {
            advMsg <- ""
            if (adversarial) advMsg <- sprintf(", disc: %.4f", as.numeric(discLoss))
            message(sprintf("Stage 1 - Epoch %d/%d - loss: %.4f, recon: %.4f, commit: %.4f%s",
                            epoch, stage1Epochs,
                            as.numeric(h$history$total_loss),
                            as.numeric(h$history$recon_loss),
                            as.numeric(h$history$commit_loss),
                            advMsg))
        }
    }

    # ========== STAGE 2: SUPERVISED VQ-VAE + CLUSTER LOSS (optional) ==========
    if (supervised) {
        stage2Epochs <- as.integer(epochs - clusterWarmup)
        message(sprintf("\n=== STAGE 2: VQ-VAE + supervised (%d epochs) ===", stage2Epochs))

        lambdaClusterFloat <- tf$constant(lambdaCluster, dtype = tf$float32)
        nClustersFloat <- tf$constant(as.numeric(nClusters), dtype = tf$float32)

        VQVAEClusterModel <- keras3::new_model_class(
            classname = "VQVAEClusterModel",
            initialize = function(enc, dec, nK, latDim, ...) {
                super$initialize(...)
                self$enc <- enc
                self$dec <- dec
                self$vq_codebook <- self$add_weight(
                    name = "vq_codebook_s2",
                    shape = as.integer(c(nK, latDim)),
                    initializer = keras3::initializer_zeros(),
                    trainable = !useEMA
                )
                if (useEMA) {
                    self$ema_count <- self$add_weight(
                        name = "ema_count_s2", shape = list(as.integer(nK)),
                        initializer = keras3::initializer_ones(), trainable = FALSE)
                    self$ema_sum <- self$add_weight(
                        name = "ema_sum_s2", shape = as.integer(c(nK, latDim)),
                        initializer = keras3::initializer_zeros(), trainable = FALSE)
                }
                self$total_loss_tracker <- keras3::metric_mean(name = "total_loss")
                self$recon_loss_tracker <- keras3::metric_mean(name = "recon_loss")
                self$cluster_loss_tracker <- keras3::metric_mean(name = "cluster_loss")
            },
            metrics = keras3::mark_active(function() {
                list(self$total_loss_tracker, self$recon_loss_tracker,
                     self$cluster_loss_tracker)
            }),
            train_step = function(data) {
                x <- data[[1]]
                yOneHot <- data[[2]]

                with(tf$GradientTape() %as% tape, {
                    z_e <- self$enc(x, training = TRUE)

                    z_e_exp <- tf$expand_dims(z_e, axis = 1L)
                    cb_exp <- tf$expand_dims(self$vq_codebook, axis = 0L)
                    distances <- tf$reduce_sum(tf$square(z_e_exp - cb_exp), axis = 2L)

                    if (softQuantize) {
                        weights <- tf$nn$softmax(-distances / temperatureFloat, axis = 1L)
                        z_q <- tf$linalg$matmul(weights, self$vq_codebook)
                    } else {
                        indices <- tf$argmin(distances, axis = 1L)
                        z_q <- tf$gather(self$vq_codebook, indices)
                    }

                    if (skipVQDecode) {
                        x_recon <- self$dec(z_e, training = TRUE)
                    } else {
                        z_q_st <- z_e + tf$stop_gradient(z_q - z_e)
                        x_recon <- self$dec(z_q_st, training = TRUE)
                    }

                    if (!is.null(reconWeightsTF)) {
                        reconLoss <- tf$reduce_mean(tf$square(x - x_recon) * reconWeightsTF) * inputDimFloat
                    } else {
                        reconLoss <- tf$reduce_mean(tf$square(x - x_recon)) * inputDimFloat
                    }
                    commitLoss <- tf$reduce_mean(tf$square(z_e - tf$stop_gradient(z_q)))

                    # Cluster loss: softmax over negative distances
                    clusterProbs <- tf$nn$softmax(-distances / temperatureFloat, axis = 1L)
                    smoothedY <- yOneHot * (1.0 - labelSmoothing) +
                        labelSmoothing / nClustersFloat
                    clusterLoss <- tf$reduce_mean(
                        tf$keras$losses$categorical_crossentropy(smoothedY, clusterProbs))

                    if (useEMA) {
                        totalLoss <- reconLoss + betaFloat * commitLoss +
                                     lambdaClusterFloat * clusterLoss
                    } else {
                        codebookLoss <- tf$reduce_mean(tf$square(tf$stop_gradient(z_e) - z_q))
                        totalLoss <- reconLoss + betaFloat * commitLoss + codebookLoss +
                                     lambdaClusterFloat * clusterLoss
                    }
                })

                grads <- tape$gradient(totalLoss, self$trainable_weights)
                self$optimizer$apply_gradients(keras3::zip_lists(grads, self$trainable_weights))

                # EMA codebook update
                if (useEMA) {
                    hardIdx <- tf$argmin(distances, axis = 1L)
                    oneHotAssign <- tf$one_hot(hardIdx, depth = nKInt)
                    n_k <- tf$reduce_sum(oneHotAssign, axis = 0L)
                    sum_k <- tf$linalg$matmul(tf$transpose(oneHotAssign), z_e)
                    self$ema_count$assign(emaDecayFloat * self$ema_count +
                                          (1.0 - emaDecayFloat) * n_k)
                    self$ema_sum$assign(emaDecayFloat * self$ema_sum +
                                        (1.0 - emaDecayFloat) * sum_k)
                    updated_cb <- self$ema_sum / tf$expand_dims(self$ema_count + 1e-5, axis = 1L)
                    self$vq_codebook$assign(updated_cb)
                }

                self$total_loss_tracker$update_state(totalLoss)
                self$recon_loss_tracker$update_state(reconLoss)
                self$cluster_loss_tracker$update_state(clusterLoss)

                list(total_loss = self$total_loss_tracker$result(),
                     recon_loss = self$recon_loss_tracker$result(),
                     cluster_loss = self$cluster_loss_tracker$result())
            }
        )

        optimizerS2 <- keras3::optimizer_adam(learning_rate = learningRate * 0.5)
        vqClusterModel <- VQVAEClusterModel(enc = encoder, dec = decoder,
                                             nK = nClusters, latDim = latentDim)
        # Transfer codebook (and EMA state) from Stage 1
        cb_stage1 <- as.matrix(vqvaeModel$vq_codebook$numpy())
        vqClusterModel$vq_codebook$assign(tf$constant(cb_stage1, dtype = tf$float32))
        if (useEMA) {
            vqClusterModel$ema_count$assign(vqvaeModel$ema_count * 1.0)
            vqClusterModel$ema_sum$assign(vqvaeModel$ema_sum * 1.0)
        }
        keras3::compile(vqClusterModel, optimizer = optimizerS2)

        history$cluster <- c()
        for (epoch in seq_len(stage2Epochs)) {
            idx <- sample(nTrain)
            xShuffled <- xTrain[idx, , drop = FALSE]
            yShuffled <- yTrainOneHot[idx, , drop = FALSE]
            xTensor <- tf$constant(xShuffled, dtype = tf$float32)
            yTensor <- tf$constant(yShuffled, dtype = tf$float32)

            dataset <- tf$data$Dataset$from_tensor_slices(
                reticulate::tuple(xTensor, yTensor))$batch(batchSize)$prefetch(tf$data$AUTOTUNE)
            h <- vqClusterModel$fit(dataset, epochs = 1L, verbose = 0L)

            history$cluster <- c(history$cluster, as.numeric(h$history$cluster_loss))

            # Dead code reset
            if (codebookReset && (epoch %% resetInterval == 0)) {
                z_all <- as.matrix(encoder(tf$constant(xTrain, dtype = tf$float32), training = FALSE))
                z_all_t <- tf$constant(z_all, dtype = tf$float32)
                z_exp <- tf$expand_dims(z_all_t, axis = 1L)
                cb_exp <- tf$expand_dims(vqClusterModel$vq_codebook, axis = 0L)
                dists <- tf$reduce_sum(tf$square(z_exp - cb_exp), axis = 2L)
                assignments <- as.integer(as.matrix(tf$argmin(dists, axis = 1L)))
                counts <- tabulate(assignments + 1L, nbins = nClusters)
                dead <- which(counts == 0)
                if (length(dead) > 0) {
                    message(sprintf("  Epoch %d: resetting %d dead codebook entries", epoch, length(dead)))
                    randIdx <- sample(nTrain, length(dead), replace = TRUE)
                    newVals <- z_all[randIdx, , drop = FALSE] +
                        matrix(rnorm(length(dead) * latentDim, sd = 0.01),
                               nrow = length(dead))
                    cbNumpy <- as.matrix(vqClusterModel$vq_codebook$numpy())
                    cbNumpy[dead, ] <- newVals
                    vqClusterModel$vq_codebook$assign(tf$constant(cbNumpy, dtype = tf$float32))
                    if (useEMA) {
                        emaSumNp <- as.matrix(vqClusterModel$ema_sum$numpy())
                        emaCntNp <- as.numeric(vqClusterModel$ema_count$numpy())
                        emaSumNp[dead, ] <- newVals
                        emaCntNp[dead] <- 1.0
                        vqClusterModel$ema_sum$assign(tf$constant(emaSumNp, dtype = tf$float32))
                        vqClusterModel$ema_count$assign(tf$constant(emaCntNp, dtype = tf$float32))
                    }
                }
            }

            # Adversarial batch removal step (Stage 2)
            if (adversarial) {
                z_all_adv <- encoder(tf$constant(xTrain[idx, , drop = FALSE], dtype = tf$float32), training = FALSE)
                sLabels <- tf$constant(sampleLabelsTrain[idx], dtype = tf$int32)
                sOneHot <- tf$one_hot(sLabels, depth = as.integer(nSamplesAdv))

                for (dStep in 1:3) {
                    with(tf$GradientTape() %as% dtape, {
                        discPreds <- discriminator(tf$stop_gradient(z_all_adv), training = TRUE)
                        discLoss <- tf$reduce_mean(
                            tf$keras$losses$categorical_crossentropy(sOneHot, discPreds))
                    })
                    dgrads <- dtape$gradient(discLoss, discriminator$trainable_weights)
                    optimizerDisc$apply_gradients(
                        keras3::zip_lists(dgrads, discriminator$trainable_weights))
                }

                with(tf$GradientTape() %as% etape, {
                    z_adv <- encoder(tf$constant(xTrain[idx, , drop = FALSE], dtype = tf$float32), training = TRUE)
                    discPreds2 <- discriminator(z_adv, training = FALSE)
                    uniformTarget <- tf$ones_like(discPreds2) / nSamplesAdvFloat
                    advLoss <- -advLambdaFloat * tf$reduce_mean(
                        tf$keras$losses$categorical_crossentropy(uniformTarget, discPreds2))
                })
                egrads <- etape$gradient(advLoss, encoder$trainable_weights)
                validGrads <- list()
                validWeights <- list()
                for (gi in seq_along(egrads)) {
                    if (!is.null(egrads[[gi]])) {
                        validGrads[[length(validGrads) + 1]] <- egrads[[gi]]
                        validWeights[[length(validWeights) + 1]] <- encoder$trainable_weights[[gi]]
                    }
                }
                if (length(validGrads) > 0) {
                    optimizerAE$apply_gradients(keras3::zip_lists(validGrads, validWeights))
                }
            }

            if (epoch %% 10 == 0 || epoch == 1) {
                advMsg <- ""
                if (adversarial) advMsg <- sprintf(", disc: %.4f", as.numeric(discLoss))
                message(sprintf("Stage 2 - Epoch %d/%d - total: %.4f, cluster: %.4f%s",
                                epoch, stage2Epochs,
                                as.numeric(h$history$total_loss),
                                as.numeric(h$history$cluster_loss),
                                advMsg))
            }
        }
    }

    message("\nTraining complete!")

    # ========== BUILD INFERENCE WRAPPERS ==========
    # Snapshot codebook from the final model
    if (supervised) {
        codebookNumpy <- as.matrix(vqClusterModel$vq_codebook$numpy())
    } else {
        codebookNumpy <- as.matrix(vqvaeModel$vq_codebook$numpy())
    }

    # ClusterDist Layer: z_e -> softmax(-distances) -> cluster probabilities
    ClusterDistLayer <- keras3::new_layer_class(
        classname = "ClusterDistLayer",
        initialize = function(cb_init, ...) {
            super$initialize(...)
            self$cb <- self$add_weight(
                name = "cluster_codebook", shape = dim(cb_init),
                initializer = keras3::initializer_constant(value = cb_init),
                trainable = FALSE)
        },
        call = function(z_e) {
            z_exp <- tf$expand_dims(z_e, axis = 1L)
            cb_exp <- tf$expand_dims(self$cb, axis = 0L)
            dists <- tf$reduce_sum(tf$square(z_exp - cb_exp), axis = 2L)
            tf$nn$softmax(-dists, axis = 1L)
        }
    )

    if (skipVQDecode) {
        # Decoder takes z_e directly (no VQ layer needed)
        wrappedDecoder <- decoder
    } else {
        # VQ Layer: z_e -> nearest codebook vector -> z_q
        VQLayer <- keras3::new_layer_class(
            classname = "VQLayer",
            initialize = function(cb_init, ...) {
                super$initialize(...)
                self$cb <- self$add_weight(
                    name = "vq_codebook", shape = dim(cb_init),
                    initializer = keras3::initializer_constant(value = cb_init),
                    trainable = FALSE)
            },
            call = function(z_e) {
                z_exp <- tf$expand_dims(z_e, axis = 1L)
                cb_exp <- tf$expand_dims(self$cb, axis = 0L)
                dists <- tf$reduce_sum(tf$square(z_exp - cb_exp), axis = 2L)
                indices <- tf$argmin(dists, axis = 1L)
                tf$gather(self$cb, indices)
            }
        )
        vqLayer <- VQLayer(cb_init = codebookNumpy)
        wrapDecInput <- keras3::layer_input(shape = latentDim, name = "vq_dec_input")
        z_q_wrap <- vqLayer(wrapDecInput)
        x_recon_wrap <- decoder(z_q_wrap)
        wrappedDecoder <- keras3::keras_model(wrapDecInput, x_recon_wrap, name = "vq_decoder")
    }

    # Cluster head: input z_e -> softmax(-distances)
    clusterDistLayer <- ClusterDistLayer(cb_init = codebookNumpy)
    clInput <- keras3::layer_input(shape = latentDim, name = "cluster_input")
    clProbs <- clusterDistLayer(clInput)
    clusterHead <- keras3::keras_model(clInput, clProbs, name = "cluster_head")

    # Log codebook utilization
    z_final <- as.matrix(encoder(tf$constant(xTrain, dtype = tf$float32), training = FALSE))
    z_final_t <- tf$constant(z_final, dtype = tf$float32)
    z_exp <- tf$expand_dims(z_final_t, axis = 1L)
    cb_final <- tf$constant(codebookNumpy, dtype = tf$float32)
    cb_exp <- tf$expand_dims(cb_final, axis = 0L)
    dists_final <- tf$reduce_sum(tf$square(z_exp - cb_exp), axis = 2L)
    final_assignments <- as.integer(as.matrix(tf$argmin(dists_final, axis = 1L)))
    final_counts <- tabulate(final_assignments + 1L, nbins = nClusters)
    message(sprintf("  Codebook utilization: %d/%d entries used", sum(final_counts > 0), nClusters))
    message(sprintf("  Cells per entry: %s", paste(final_counts, collapse = ", ")))

    result <- list(
        encoder = encoder,
        decoder = wrappedDecoder,
        clusterHead = clusterHead,
        rawDecoder = decoder,
        codebook = codebookNumpy,
        history = history,
        nClusters = nClusters,
        latentDim = latentDim,
        supervised = supervised
    )

    return(result)
}


#' Compute Cell Importance Scores
#'
#' Computes importance weights for each cell based on reconstruction quality.
#' Cells with lower reconstruction error are considered more "representative"
#' and receive higher importance weights. This can be used for attention-weighted
#' feature aggregation.
#'
#' @param data Matrix of marker expression data (cells x markers).
#' @param encoder Trained encoder model from trainSWAEClustering.
#' @param decoder Trained decoder model from trainSWAEClustering.
#' @param method Method for computing importance: "recon_error" (inverse reconstruction error),
#'   "cluster_confidence" (max cluster probability), or "combined" (product of both).
#' @param clusterHead Optional cluster head model for cluster_confidence method.
#' @param temperature Temperature for softmax conversion of scores. Lower = more peaked.
#' @param normalize If TRUE, normalize weights to sum to 1 within each sample.
#' @return Numeric vector of importance weights for each cell.
#' @export
#' @examples
#' \dontrun{
#' model <- trainSWAEClustering(trainData, useMarkers, nClusters = 10)
#' importance <- computeCellImportance(
#'   data = allData[, useMarkers],
#'   encoder = model$encoder,
#'   decoder = model$decoder,
#'   method = "recon_error"
#' )
#' }
computeCellImportance <- function(data, encoder, decoder, method = "recon_error",
                                   clusterHead = NULL, temperature = 1.0,
                                   normalize = FALSE) {
    data <- as.matrix(data)
    n <- nrow(data)

    # Get latent representations and reconstructions
    latent <- as.matrix(encoder(data))
    recon <- as.matrix(decoder(latent))

    if (method == "recon_error" || method == "combined") {
        reconError <- rowSums((data - recon)^2)
        reconScores <- -reconError
        reconScores <- reconScores - max(reconScores)
        reconWeights <- exp(reconScores / temperature)
    }

    if (method == "cluster_confidence" || method == "combined") {
        if (is.null(clusterHead)) {
            stop("clusterHead required for cluster_confidence method")
        }
        clusterProbs <- as.matrix(clusterHead(latent))
        confScores <- apply(clusterProbs, 1, max)
        confWeights <- confScores
    }

    if (method == "recon_error") {
        weights <- reconWeights
    } else if (method == "cluster_confidence") {
        weights <- confWeights
    } else if (method == "combined") {
        weights <- reconWeights * confWeights
    } else {
        stop("Unknown method. Use 'recon_error', 'cluster_confidence', or 'combined'")
    }

    if (normalize) {
        weights <- weights / sum(weights)
    }

    return(weights)
}


#' Unsupervised Contrastive Loss (SimCLR-style)
#'
#' Computes contrastive loss using augmented views of the same cell.
#' Pulls together augmented versions of the same cell, pushes apart different cells.
#'
#' @param z1 Tensor of embeddings for first augmented view.
#' @param z2 Tensor of embeddings for second augmented view.
#' @param temperature Temperature parameter for softmax.
#' @return Scalar contrastive loss.
#' @noRd
unsupervisedContrastiveLoss <- function(z1, z2, temperature = 0.5) {
    tf <- tensorflow::tf

    z1Norm <- tf$math$l2_normalize(z1, axis = 1L)
    z2Norm <- tf$math$l2_normalize(z2, axis = 1L)

    batchSize <- tf$shape(z1)[1]
    zAll <- tf$concat(list(z1Norm, z2Norm), axis = 0L)

    similarity <- tf$matmul(zAll, zAll, transpose_b = TRUE) / temperature

    labels <- tf$concat(list(
        tf$range(batchSize, 2L * batchSize),
        tf$range(0L, batchSize)
    ), axis = 0L)

    mask <- tf$ones(c(2L * batchSize, 2L * batchSize)) - tf$eye(2L * batchSize)
    similarity <- similarity * mask - (1 - mask) * 1e9

    loss <- tf$reduce_mean(
        tf$nn$sparse_softmax_cross_entropy_with_logits(labels = labels, logits = similarity)
    )

    return(loss)
}


#' Data Augmentation for Cells
#'
#' Applies random augmentation to cell marker data for contrastive learning.
#' Augmentations include: Gaussian noise, dropout (masking), and scaling.
#'
#' @param x Tensor of cell data.
#' @param noiseScale Scale of Gaussian noise to add.
#' @param dropoutRate Probability of masking each feature.
#' @param scaleRange Range for random scaling factor.
#' @return Augmented tensor.
#' @noRd
augmentCells <- function(x, noiseScale = 0.1, dropoutRate = 0.1, scaleRange = c(0.9, 1.1)) {
    tf <- tensorflow::tf

    noise <- tf$random$normal(tf$shape(x), mean = 0.0, stddev = noiseScale)
    xAug <- x + noise

    mask <- tf$cast(tf$random$uniform(tf$shape(x)) > dropoutRate, tf$float32)
    xAug <- xAug * mask

    scale <- tf$random$uniform(shape = list(), minval = scaleRange[1], maxval = scaleRange[2])
    xAug <- xAug * scale

    xAug <- tf$clip_by_value(xAug, 0.0, 1.0)

    return(xAug)
}


#' Compute Sliced Wasserstein Barycenter
#'
#' Computes the Wasserstein barycenter of multiple samples using the sliced
#' Wasserstein distance. The barycenter is a representative distribution that
#' minimizes the average transport cost to all input samples. This is used as
#' the reference distribution for OT batch correction.
#'
#' @param samplesList List of data frames or matrices, one per sample.
#'   Each element contains cells (rows) x markers (columns).
#' @param nCellsPerSample Integer; number of cells to subsample from each sample.
#'   Default is 2000.
#' @param nBarycenterCells Integer; number of cells in the output barycenter.
#'   Default is 5000.
#' @param nProjections Integer; number of random projections per iteration.
#'   Default is 100.
#' @param nIter Integer; number of barycenter iterations. Default is 10.
#' @param seed Integer; random seed. Default is 1994.
#' @return Matrix of barycenter cells (nBarycenterCells x markers).
#' @export
#' @examples
#' \dontrun{
#' # Split decoded data by sample
#' decoded_df <- as.data.frame(decoded_data)
#' decoded_df$sample_id <- df$sample_id
#' samples_list <- split(decoded_df[, markers], decoded_df$sample_id)
#'
#' # Compute barycenter
#' barycenter <- computeSlicedBarycenter(
#'   samplesList = samples_list,
#'   nCellsPerSample = 10000,
#'   nBarycenterCells = 5000,
#'   nProjections = 200,
#'   nIter = 20
#' )
#' }
computeSlicedBarycenter <- function(samplesList, nCellsPerSample = 2000L,
                                     nBarycenterCells = 5000L,
                                     nProjections = 100L, nIter = 10L,
                                     seed = 1994) {
    set.seed(seed)

    nSamples <- length(samplesList)
    nCellsPerSample <- as.integer(nCellsPerSample)
    nBarycenterCells <- as.integer(nBarycenterCells)
    nProjections <- as.integer(nProjections)
    nIter <- as.integer(nIter)

    message(sprintf("Computing sliced Wasserstein barycenter:"))
    message(sprintf("  %d samples, %d cells per sample, %d barycenter cells",
                   nSamples, nCellsPerSample, nBarycenterCells))
    message(sprintf("  %d projections, %d iterations", nProjections, nIter))

    # Subsample each distribution to equal size
    subsampledList <- lapply(samplesList, function(x) {
        x <- as.matrix(x)
        n <- nrow(x)
        idx <- sample(n, min(nCellsPerSample, n), replace = nCellsPerSample > n)
        x[idx, , drop = FALSE]
    })

    d <- ncol(subsampledList[[1]])  # marker dimension

    # Initialize barycenter from pooled subsampled cells
    allCells <- do.call(rbind, subsampledList)
    baryIdx <- sample(nrow(allCells), nBarycenterCells,
                      replace = nBarycenterCells > nrow(allCells))
    barycenter <- allCells[baryIdx, , drop = FALSE]

    # Iterative barycenter computation via sliced OT (Bonneel et al. 2015)
    # Key idea: work entirely in 1D per projection, then back-project
    # displacements along the projection direction
    for (iter in seq_len(nIter)) {
        # Generate random projection directions (unit vectors)
        theta <- matrix(stats::rnorm(d * nProjections), nrow = d, ncol = nProjections)
        theta <- sweep(theta, 2, sqrt(colSums(theta^2)), "/")

        # Accumulate displacement in d-dimensional space
        displacement <- matrix(0, nrow = nBarycenterCells, ncol = d)

        for (p in seq_len(nProjections)) {
            direction <- theta[, p]

            # Project barycenter onto 1D
            baryProj <- as.vector(barycenter %*% direction)
            baryOrder <- order(baryProj)
            barySorted1D <- baryProj[baryOrder]

            # Compute average sorted 1D projection across all samples
            avg1D <- numeric(nBarycenterCells)

            for (s in seq_len(nSamples)) {
                sampleData <- subsampledList[[s]]
                sampleProj <- as.vector(sampleData %*% direction)
                sampleSorted1D <- sort(sampleProj)
                nSampleCells <- length(sampleSorted1D)

                # Interpolate sample's sorted 1D projections to barycenter size
                sourceQuantiles <- seq(0, 1, length.out = nSampleCells)
                targetQuantiles <- seq(0, 1, length.out = nBarycenterCells)
                matched1D <- stats::approx(sourceQuantiles, sampleSorted1D,
                                           xout = targetQuantiles, rule = 2)$y
                avg1D <- avg1D + matched1D
            }

            avg1D <- avg1D / nSamples

            # 1D displacement: where the barycenter points should move to
            delta1D <- avg1D - barySorted1D

            # Back-project 1D displacement to d dimensions along direction
            # displacement[baryOrder[i], ] += delta1D[i] * direction
            displacement[baryOrder, ] <- displacement[baryOrder, ] +
                tcrossprod(delta1D, direction)
        }

        # Average displacement over projections and update barycenter
        barycenter <- barycenter + displacement / nProjections

        if (iter %% max(1, nIter %/% 5) == 0 || iter == 1) {
            message(sprintf("  Iteration %d/%d", iter, nIter))
        }
    }

    # Clip to valid range
    barycenter <- pmax(barycenter, 0)
    barycenter <- pmin(barycenter, 1)

    colnames(barycenter) <- colnames(subsampledList[[1]])
    message("Barycenter computation complete")

    return(barycenter)
}


#' Compute Sliced Wasserstein Barycenter via Autograd
#'
#' Computes the barycenter (Frechet mean) of multiple empirical distributions
#' using gradient descent on the sliced Wasserstein-2 distance, with TensorFlow
#' automatic differentiation. This is a more principled alternative to the
#' fixed-point iteration in \code{\link{computeSlicedBarycenter}}.
#'
#' @param samplesList List of data frames or matrices, one per sample. Each
#'   element has cells as rows and markers as columns.
#' @param nCellsPerSample Integer; number of cells to subsample from each
#'   distribution. Default is 2000.
#' @param nBarycenterCells Integer; number of cells in the barycenter. Default
#'   is 5000.
#' @param nProjections Integer; number of random projections per gradient step.
#'   Default is 50.
#' @param nIter Integer; number of gradient descent iterations. Default is 200.
#' @param learningRate Numeric; Adam learning rate. Default is 0.01.
#' @param sampleBatchSize Integer; number of samples per mini-batch. Processes
#'   samples in chunks to avoid OOM. Default is 16.
#' @param seed Integer; random seed. Default is 1994.
#' @return Matrix of barycenter cells (nBarycenterCells x markers).
#' @export
#' @examples
#' \dontrun{
#' barycenter <- computeSlicedBarycenterAutograd(
#'   samplesList = samples_list,
#'   nCellsPerSample = 10000,
#'   nBarycenterCells = 5000,
#'   nProjections = 50,
#'   nIter = 200
#' )
#' }
computeSlicedBarycenterAutograd <- function(samplesList, nCellsPerSample = 2000L,
                                             nBarycenterCells = 5000L,
                                             nProjections = 50L, nIter = 200L,
                                             learningRate = 0.01,
                                             sampleBatchSize = 16L,
                                             seed = 1994) {
    set.seed(seed)
    tf$random$set_seed(as.integer(seed))

    nSamples <- length(samplesList)
    nBarycenterCells <- as.integer(nBarycenterCells)
    nProjections <- as.integer(nProjections)
    nIter <- as.integer(nIter)
    sampleBatchSize <- as.integer(min(sampleBatchSize, nSamples))

    message(sprintf("Computing sliced Wasserstein barycenter (autograd):"))
    message(sprintf("  %d samples, %d cells per sample, %d barycenter cells",
                   nSamples, nCellsPerSample, nBarycenterCells))
    message(sprintf("  %d projections/step, %d iterations, lr=%.4f, batch=%d samples",
                   nProjections, nIter, learningRate, sampleBatchSize))

    # Subsample each distribution to exactly nBarycenterCells (same size for sort)
    subsampledList <- lapply(samplesList, function(x) {
        x <- as.matrix(x)
        n <- nrow(x)
        idx <- sample(n, nBarycenterCells, replace = nBarycenterCells > n)
        x[idx, , drop = FALSE]
    })

    d <- ncol(subsampledList[[1]])

    # Convert each sample to a TF constant (keep as list to avoid huge tensor)
    samplesTF <- unname(lapply(subsampledList, function(x) {
        tf$constant(x, dtype = tf$float32)
    }))

    # Initialize barycenter from pooled subsampled cells
    allCells <- do.call(rbind, subsampledList)
    baryIdx <- sample(nrow(allCells), nBarycenterCells,
                      replace = nBarycenterCells > nrow(allCells))
    baryInit <- allCells[baryIdx, , drop = FALSE]

    barycenter <- tf$Variable(baryInit, dtype = tf$float32)

    # Adam optimizer
    optimizer <- tf$keras$optimizers$Adam(learning_rate = learningRate)

    for (iter in seq_len(nIter)) {
        # Random projections on unit sphere: (d, nProjections)
        theta <- tf$random$normal(shape = c(as.integer(d), nProjections))
        theta <- theta / tf$norm(theta, axis = 0L, keepdims = TRUE)

        # Mini-batch of samples (stochastic over samples)
        batchIdx <- sample(nSamples, sampleBatchSize)

        # Stack mini-batch into tensor: (B, nBarycenterCells, d)
        batchTensor <- tf$stack(samplesTF[batchIdx])

        with(tf$GradientTape() %as% tape, {
            # Project barycenter: (nBarycenterCells, L)
            baryProj <- tf$matmul(barycenter, theta)
            # Expand for broadcasting: (1, nBarycenterCells, L)
            baryProj <- tf$expand_dims(baryProj, 0L)

            # Project mini-batch: (B, nBarycenterCells, L)
            batchProj <- tf$einsum("snd,dL->snL", batchTensor, theta)

            # Sort along cell axis (axis=1)
            barySorted <- tf$sort(baryProj, axis = 1L)
            batchSorted <- tf$sort(batchProj, axis = 1L)

            # SW_2 distance: average over batch, cells, projections
            sw2 <- tf$reduce_mean(tf$square(barySorted - batchSorted))
        })

        grads <- tape$gradient(sw2, barycenter)
        optimizer$apply_gradients(list(list(grads, barycenter)))

        # Clip to [0, 1]
        barycenter$assign(tf$clip_by_value(barycenter, 0.0, 1.0))

        if (iter %% max(1L, nIter %/% 5L) == 0L || iter == 1L) {
            loss_val <- as.numeric(sw2)
            message(sprintf("  Iteration %d/%d - SW2 loss: %.6f",
                           iter, nIter, loss_val))
        }
    }

    result <- as.matrix(barycenter$numpy())
    colnames(result) <- colnames(subsampledList[[1]])
    message("Barycenter computation complete (autograd)")

    return(result)
}


#' Batch Correct Using Sliced Optimal Transport
#'
#' Corrects batch effects by transporting a target distribution toward a reference
#' distribution using sliced optimal transport. The interpolation parameter alpha
#' controls the strength of correction (0 = no correction, 1 = full transport).
#'
#' @param targetData Matrix of target data to correct (cells x markers).
#' @param reference Matrix of reference distribution (e.g., barycenter from
#'   computeSlicedBarycenter).
#' @param alpha Numeric; interpolation parameter between 0 (no correction) and
#'   1 (full transport to reference). Default is 0.7.
#' @param nProjections Integer; number of random projections. Default is 100.
#' @param seed Integer; random seed. Default is 1994.
#' @return List containing:
#'   \item{corrected}{Matrix of corrected data (same dimensions as targetData).}
#'   \item{alpha}{The alpha value used.}
#' @export
#' @examples
#' \dontrun{
#' # Correct target data toward the barycenter
#' result <- batchCorrectOT(
#'   targetData = decoded_data[, markers],
#'   reference = barycenter,
#'   alpha = 0.7,
#'   nProjections = 200
#' )
#' corrected <- result$corrected
#' }
batchCorrectOT <- function(targetData, reference, alpha = 0.7, nProjections = 100L,
                            seed = 1994) {
    set.seed(seed)

    targetData <- as.matrix(targetData)
    reference <- as.matrix(reference)

    nTarget <- nrow(targetData)
    nRef <- nrow(reference)
    d <- ncol(targetData)

    nProjections <- as.integer(nProjections)

    message(sprintf("OT batch correction: %d target cells, %d reference cells, alpha=%.2f",
                   nTarget, nRef, alpha))

    # Generate random projection directions
    theta <- matrix(stats::rnorm(d * nProjections), nrow = d, ncol = nProjections)
    theta <- sweep(theta, 2, sqrt(colSums(theta^2)), "/")

    # Accumulate displacement across all projections
    displacement <- matrix(0, nrow = nTarget, ncol = d)

    for (p in seq_len(nProjections)) {
        direction <- theta[, p]

        # Project target and reference onto 1D
        targetProj <- as.vector(targetData %*% direction)
        refProj <- as.vector(reference %*% direction)

        # Sort both projections
        targetOrder <- order(targetProj)
        targetProjSorted <- targetProj[targetOrder]
        refSorted <- sort(refProj)

        # 1D OT map: T(x) = F_ref^{-1}(F_target(x))
        # Interpolate reference sorted values to match target quantiles
        targetQuantiles <- seq(0, 1, length.out = nTarget)
        refQuantiles <- seq(0, 1, length.out = nRef)
        refInterp <- stats::approx(refQuantiles, refSorted,
                                    xout = targetQuantiles, rule = 2)$y

        # 1D displacement for sorted target cells
        delta1D <- refInterp - targetProjSorted

        # Back-project 1D displacement to d dimensions along direction
        displacement[targetOrder, ] <- displacement[targetOrder, ] +
            tcrossprod(delta1D, direction)
    }

    # Average over projections, apply alpha interpolation
    corrected <- targetData + alpha * displacement / nProjections

    # Clip to valid range
    corrected <- pmax(corrected, 0)
    corrected <- pmin(corrected, 1)

    colnames(corrected) <- colnames(targetData)
    message("OT batch correction complete")

    return(list(corrected = corrected, alpha = alpha))
}


#' Per-Cluster OT Batch Correction
#'
#' Corrects batch effects by applying sliced OT correction independently to each
#' cell cluster. Computes a cluster-specific barycenter from per-sample splits,
#' then transports each cluster's cells to its own reference distribution.
#'
#' @param decodedData Matrix of decoded marker expression (cells x markers).
#' @param sampleIds Character vector of sample IDs for each cell.
#' @param clusters Integer vector of cluster assignments for each cell.
#' @param alpha Interpolation strength (0 = no correction, 1 = full transport).
#' @param nProjections Number of random projections for sliced OT.
#' @param nBarycenterCells Number of cells in each cluster-specific barycenter.
#' @param nBarycenterIter Number of autograd barycenter iterations.
#' @param minCellsPerCluster Minimum cells per cluster per sample to include in barycenter.
#' @param seed Random seed.
#' @return List with \code{corrected} matrix (same dims as input).
#' @export
batchCorrectOTPerCluster <- function(decodedData, sampleIds, clusters,
                                      alpha = 0.7, nProjections = 200L,
                                      nBarycenterCells = 5000L,
                                      nBarycenterIter = 100L,
                                      minCellsPerCluster = 50L,
                                      seed = 1994) {
    set.seed(seed)
    decodedData <- as.matrix(decodedData)
    markers <- colnames(decodedData)
    corrected <- decodedData  # copy
    uniqueClusters <- sort(unique(clusters))

    message(sprintf("Per-cluster OT correction: %d clusters, %d cells, alpha=%.2f",
                    length(uniqueClusters), nrow(decodedData), alpha))

    for (k in uniqueClusters) {
        clusterIdx <- which(clusters == k)
        nClusterCells <- length(clusterIdx)

        if (nClusterCells < 100) {
            message(sprintf("  Cluster %s: %d cells (skipped, too few)", k, nClusterCells))
            next
        }

        clusterData <- decodedData[clusterIdx, , drop = FALSE]
        clusterSamples <- sampleIds[clusterIdx]

        # Build per-sample splits for this cluster
        clusterDf <- as.data.frame(clusterData)
        colnames(clusterDf) <- markers
        clusterDf$sample_id <- clusterSamples
        samplesList <- split(clusterDf[, markers], clusterDf$sample_id)

        # Filter out samples with too few cells in this cluster
        samplesList <- samplesList[sapply(samplesList, nrow) >= minCellsPerCluster]

        if (length(samplesList) < 3) {
            message(sprintf("  Cluster %s: %d cells, <3 samples with enough cells (skipped)",
                            k, nClusterCells))
            next
        }

        # Compute cluster-specific barycenter
        nBary <- min(nBarycenterCells, min(sapply(samplesList, nrow)))
        barycenter <- computeSlicedBarycenterAutograd(
            samplesList = samplesList,
            nCellsPerSample = nBary,
            nBarycenterCells = nBary,
            nProjections = 50L,
            nIter = nBarycenterIter,
            learningRate = 0.01,
            sampleBatchSize = min(16L, length(samplesList)),
            seed = seed
        )

        # OT correct this cluster
        result <- batchCorrectOT(
            targetData = clusterData,
            reference = barycenter,
            alpha = alpha,
            nProjections = nProjections,
            seed = seed
        )
        corrected[clusterIdx, ] <- result$corrected

        message(sprintf("  Cluster %s: %d cells corrected", k, nClusterCells))
    }

    colnames(corrected) <- markers
    message("Per-cluster OT correction complete")
    return(list(corrected = corrected, alpha = alpha))
}
