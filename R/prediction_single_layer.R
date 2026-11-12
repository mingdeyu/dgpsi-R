# The Python linked predictor assumes a first layer followed by a final layer.
# Evaluate a one-layer network directly through its native component predictors.
# Moment aggregation includes both conditional variance and between-imputation
# variation, matching the linked predictor's mixture calculation.
.predict_single_layer_lgp <- function(engine, x, method, full_layer, sample_size, M) {
  helper <- reticulate::py_run_string("
import numpy as np
def single_layer_predict(engine, x, method, full_layer, sample_size, m):
    import numpy as np
    means, variances, samples = [], [], []
    for layers in engine.all_layer_set:
        mm, vv, ss = [], [], []
        for model in layers[0]:
            inputs = x[:, model.local_input_idx]
            if model.type == 'gp':
                mu, var = engine.gp_pred(inputs, None, None, None, model.structure, m)
            else:
                _, _, mu, var = engine.dgp_pred(inputs, None, None, None, model.structure, m)
            mm.append(mu)
            vv.append(var)
            if method == 'sampling':
                draw = np.random.normal(mu, np.sqrt(np.maximum(var, 0)),
                                        size=(sample_size,) + mu.shape)
                ss.append(draw.transpose(2, 1, 0))
        means.append(mm)
        variances.append(vv)
        samples.append(ss)
    if method == 'sampling':
        result = [np.concatenate([s[k] for s in samples], axis=2)
                  for k in range(len(samples[0]))]
        return [result] if full_layer else result
    mu = [np.mean([a[k] for a in means], axis=0) for k in range(len(means[0]))]
    var = [np.maximum(np.mean([v[k] + a[k]**2 for a, v in zip(means, variances)],
                              axis=0) - mu[k]**2, 0) for k in range(len(mu))]
    return [[mu], [var]] if full_layer else [mu, var]
", local = TRUE)
  helper$single_layer_predict(engine, x, method, full_layer,
                              as.integer(sample_size), as.integer(M))
}
