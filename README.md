# Visual Search with heterogenous distractors
This repository accompanies the following paper, currently under review:
 [Mihali, A, Ma, WJ. The psychophysics of visual search with heterogeneous distractors, 2020, bioRxiv](https://www.biorxiv.org/content/10.1101/2020.08.10.244707v3)


## Data

We concatenated all experimental sessions (4) and blocks (8/session) into `alldata_exp1.mat` and respectively `alldata_exp2.mat`, which are Nsubj * Ncond structs with 1 field. There are 4 conditions: 

- Perception Detection
- Perception Localization
- Memory Detection
- Memory Localization


All analyses were performed on the `alldata` structs. Within each subject and condition, there is the `data` struct with 8 fields, each with variables of length Ntrials (800):


- setsize N: in [2,3,4,6]
- stims: in (-pi, pi), value of the stimuli array
- target_pres: 1 or 0 for Det, always 1 for Loc
- target_loc: 1-6 or Nan for Det, 1-6 for Loc. clockwise direction starting with loc 1 as vertical
- target_val: in (-pi, pi), value of the foveally presented target
- response: Present/Absent for Det, a location in 1-6 for Loc
- reaction_time (sec)
- sess_and_block_num: session 1-4, block 1-8


## Scripts


- We fit the optimal model both by itself and with decision noise (model 1 and model 2) for both experiments and all conditions via `TD_TL_model_fitting.m`, which finds the parameters that maximize that probability of the data (likelihood) via `Loglike_all.m`. We call  `TD_TL_model_fitting.m` on a cluster computer via `TD_TL_model_fitting.sh`.
- As the loglikelihood is a noisy function with 5-8 parameters, we made use of an optimizer that is better suited for navigating noisy landscapes [Bayesian adaptive direct search (BADS) algorithm](https://github.com/lacerbi/bads). The optimization algorithm works better with several starting points to ensure that we indeed managed to find the best fitting model parameters that maximize the loglikelihood. Thus, the results from the previous script represent 20 runs for each fitting, from which we will pick the parameters that yield the highest loglikelihood. The script `fitting_pick_best.m` accomplishes this. The negative loglikelihood values and the best fitting parameters are stored in `nll_params_best_all.mat`.
- Knowing the parameters that maximize the loglikelihood for each fitting, we used them to generate the predictions of each model via `TD_TL_model_pred.m`, which in turn calls `Predict_all.m`, these being run on the cluster computer via `TD_TL_model_pred.sh`. 
- `analysis_and_plots.m` outputs the summary statistics and psychometric curves from `alldata.mat`. If the flag `model_pred` is turned on, the script loads the model predictions for each model and each subject. 



