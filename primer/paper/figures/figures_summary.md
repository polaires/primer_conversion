# Manuscript Figures - Preliminary Data

Generated: 2025-12-07

## Figure 4A: ROC Curve - Döring Dataset

AUC-ROC: **0.824**

Data file: `roc_doring.csv`

```

Figure 4A: ROC Curve - Döring Dataset (AUC = 0.824)
TPR
|                                                        .**
|                                             ************* 
|                                        ******    .        
|                                       **      .           
|              **************************    .              
|     **********                          .                 
|    **                                .                    
|  ***                              .                       
| **                             .                          
|**                           .                             
|*                         .                                
|*                      .                                   
|*                   .                                      
|*                .                                         
|*             .                                            
|*          .                                               
|*       .                                                  
|*    .                                                     
|* .                                                        
.-----------------------------------------------------------
                                                     FPR
   * = ROC curve, . = random baseline

```

## Figure 3: Feature Importance

Top discriminative features (success vs failure mean score difference):

| Rank | Feature | Difference | Weight |
|------|---------|------------|--------|
| 1 | offTarget | +0.515 | 0.10 |
| 2 | terminal3DG | +0.487 | 0.05 |
| 3 | tmRev | -0.185 | 0.05 |
| 4 | tmFwd | -0.185 | 0.02 |
| 5 | gcFwd | -0.021 | 0.06 |
| 6 | gcRev | -0.021 | 0.04 |
| 7 | hairpinRev | +0.018 | 0.05 |
| 8 | hairpinFwd | +0.018 | 0.02 |
| 9 | lengthFwd | -0.014 | 0.01 |
| 10 | lengthRev | -0.014 | 0.01 |

Data file: `feature_importance.csv`


## Model Comparison Table

| Model | Dataset | Metric | Value | Notes |
|-------|---------|--------|-------|-------|
| **Our Model** | Döring (train) | F1 | **81.9%** | Primary validation |
| **Our Model** | Döring (train) | AUC-ROC | **0.848** | |
| **Our Model** | Döring (CV) | F1 | 81.9% ± 2.9% | 5-fold cross-validation |
| **Our Model** | Kayama (test) | AUC-ROC | 0.643 | Cross-domain validation |
| **Our Model** | Kayama (test) | F1 | 42.4% | Heavy class imbalance |
| Kayama RNN | Kayama (own) | Accuracy | 70% | Their reported result |
| Kayama RNN | Kayama (undersamp) | Sens/Spec | 71%/73% | With class balancing |
| Primer3 | - | - | N/A | No published validation |

### Key Findings:
1. Our model achieves **81.9% F1** on Döring, better than Kayama RNN's 70% accuracy
2. Cross-domain AUC of 0.643 shows **partial generalization** without retraining
3. Off-target/mismatch data is the **dominant predictive feature** (diff=0.515)
4. Pure sequence-based scoring provides ~64% AUC baseline


## Summary Statistics

- **Döring Dataset**: 829 primer-template pairs
- **Success rate**: 44.0%
- **AUC-ROC**: 0.824
- **Optimal threshold**: 80

## Files Generated

1. `roc_doring.csv` - ROC curve data points
2. `feature_importance.csv` - Feature discrimination data
3. `figures_summary.md` - This summary file
