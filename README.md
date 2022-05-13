# Nonlinear CLS

Although the common linear <img src="https://render.githubusercontent.com/render/math?math=k-\epsilon"> model is utilized widely in commercial and non-commercial codes but still have some weakness in the prediction of the anisotropy of turbulence. In another hand, The cost of computation for higher-order models has made them unuseful for complex and industrial problems. Cost and benefits showed that non-linear models such as [Craft-Launder-Suga (CLS)](https://www.sciencedirect.com/science/article/abs/pii/0142727X95000796) are an appropriate choice. Here an in-house code is written for linear and non-linear CLS <img src="https://render.githubusercontent.com/render/math?math=k-\epsilon"> models to solve channel flow in Re=5600 and 14000. The results below show the performance of the written code.

| Linear model | Nonlinear model |
| --- | --- |
| <img src="https://github.com/Vaezi92/Nonlinear-CLS/blob/main/Figs/linearuplus.png" width="400"> | <img src="https://github.com/Vaezi92/Nonlinear-CLS/blob/main/Figs/nonlinearuplus.png" width="400"> |
| <img src="https://github.com/Vaezi92/Nonlinear-CLS/blob/main/Figs/linearkplus.png" width="400"> | <img src="https://github.com/Vaezi92/Nonlinear-CLS/blob/main/Figs/nonlinearkplus.png" width="400"> |

## Citation
@article{craft1996development,
  title={Development and application of a cubic eddy-viscosity model of turbulence},
  author={Craft, TJ and Launder, BE and Suga, K},
  journal={International Journal of Heat and Fluid Flow},
  volume={17},
  number={2},
  pages={108--115},
  year={1996},
  publisher={Elsevier}
}
