using GaussianProcesses: GP, MatF64, predict

function cliff_face(gpA::GP, gpB::GP, sentinels::MatF64)
    pred_A = predict(gpA, sentinels; full_cov=true)
    pred_B = predict(gpB, sentinels; full_cov=true)
    μposterior = pred_B[1].-pred_A[1]
    Σposterior = pred_B[2]+pred_A[2]
    return μposterior, Σposterior
end
