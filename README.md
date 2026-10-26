This repository contains the code to implement the empirical application described in ["Estimating Counterfactual Matrix Means with Short Panel Data"](https://arxiv.org/abs/2312.07520) by Lihua Lei and Brad Ross

The shared random seed for firm clustering and bootstrap draws is `APM_SEED`,
which defaults to `9001`. To use another seed, set it before initializing:

```bash
export APM_SEED=1234
source code/init_server.sh
```

`code/env_config.R` also defaults to `9001` when sourced without `APM_SEED` in
the environment. Seeds must be integers from `0` through `2147483647`.
Individual calls to `get_default_bootstrap()` can still override `seed`.
