# 控制函数法与2SLS系数的联系

我们知道，在$\text{2SLS}$中


$$
\hat\beta_{\text{2SLS}} = (\hat{\mathrm{X}}^{\mathrm{T}}\hat{\mathrm{X}})^{-1}\hat{\mathrm{X}}^{\mathrm{T}}\mathrm{y}
$$

其中
$$
\hat{\mathrm{X}} = (\mathrm{Z}^{\mathrm{T}}\mathrm{Z})^{-1}\mathrm{Z}^{\mathrm{T}}\mathrm{X}
$$
$\mathrm{Z}$是$N \times L$的工具变量矩阵，$\mathrm{X}$是$N \times K$的原始变量矩阵，一般有$L \gt K$。

而控制函数（$\text{Control Function, CF}$）法的思路与$\text{2SLS}$相反：$\text{2SLS}$将$\mathrm{X}$对$\mathrm{Z}$回归的预测值（投影）当做”干净的“变量重新与$\mathrm{y}$做回归；而$\text{CF}$将$\mathrm{X}$对$\mathrm{Z}$回归的残差当做内生性的来源，与$\mathrm{X}$一并与$\mathrm{y}$做回归，其$\text{DGP}$为：
$$
\mathrm{y} = \begin{bmatrix}
\mathrm{X} & \mathrm{M_Z}\mathrm{X}
\end{bmatrix}
\beta_{\text{CF}} + \epsilon_{\text{CF}}
$$
其中$\mathrm{M_Z} \equiv \mathrm{I} - \mathrm{Z}(\mathrm{Z}^{\mathrm{T}}\mathrm{Z})^{-1}\mathrm{Z}^{\mathrm{T}}$，是$\mathrm{Z}$的残差矩阵。

我们也可以同时将$\beta_{\text{CF}}$分块，即：
$$
\mathrm{y} = \begin{bmatrix}
\mathrm{X} & \mathrm{M_Z}\mathrm{X}
\end{bmatrix}
\begin{bmatrix}
\beta_{\text{CF1}} \\ \beta_{\text{CF2}}
\end{bmatrix} + \epsilon_{\text{CF}}
$$
注意：虽然$\hat\beta_{\text{2SLS}}$可以通过$\mathrm{y} \sim \hat{\mathrm{X}}$的最小二乘求出，但这只是数学上的相等（$\hat{\mathrm{X}}^{\mathrm{T}}\hat{\mathrm{X}} = \hat{\mathrm{X}}^{\mathrm{T}}\mathrm{X}$），其$\text{DGP}$依然是$\mathrm{y} = \mathrm{X}
\beta_{\text{2SLS}} + \epsilon_{\text{2SLS}}$，而不是$\mathrm{y} = \hat{\mathrm{X}}
\beta_{\text{2SLS}} + \epsilon_{\text{2SLS}}$。

有趣的是，$\hat\beta_{\text{2SLS}} = \hat\beta_{\text{CF1}}$，接下来我们就要证明这一等式。

由$\text{OLS}$知
$$
\begin{aligned}
\begin{bmatrix}
\hat\beta_{\text{CF1}} \\ \hat\beta_{\text{CF2}}
\end{bmatrix} &= \begin{bmatrix}
\mathrm{X} & \mathrm{M_Z}\mathrm{X}
\end{bmatrix}^{\mathrm{T}}\begin{bmatrix}
\mathrm{X} & \mathrm{M_Z}\mathrm{X}
\end{bmatrix} \begin{bmatrix}
\mathrm{X} & \mathrm{M_Z}\mathrm{X}
\end{bmatrix}^{\mathrm{T}} \mathrm{y} \\
&= \begin{bmatrix}
\mathrm{X}^{\mathrm{T}}\mathrm{X} & \mathrm{X}^{\mathrm{T}}\mathrm{M_Z}\mathrm{X} \\
\mathrm{X}^{\mathrm{T}}\mathrm{M_Z}\mathrm{X} & \mathrm{X}^{\mathrm{T}}\mathrm{M_Z}\mathrm{X}
\end{bmatrix}^{-1}
\begin{bmatrix}
\mathrm{X}^{\mathrm{T}} \\
\mathrm{X}^{\mathrm{T}}\mathrm{M_Z} 
\end{bmatrix} \mathrm{y}
\end{aligned}
$$
这里用到了残差矩阵$\mathrm{M_Z}$的对称幂等性。

由分块逆矩阵公式知
$$
\begin{bmatrix}
\mathrm{X}^{\mathrm{T}}\mathrm{X} & \mathrm{X}^{\mathrm{T}}\mathrm{M_Z}\mathrm{X} \\
\mathrm{X}^{\mathrm{T}}\mathrm{M_Z}\mathrm{X} & \mathrm{X}^{\mathrm{T}}\mathrm{M_Z}\mathrm{X}
\end{bmatrix}^{-1} = 
\begin{bmatrix}
(\mathrm{X}^{\mathrm{T}}(\mathrm{I} - \mathrm{M_Z})\mathrm{X})^{-1} & -(\mathrm{X}^{\mathrm{T}}(\mathrm{I} - \mathrm{M_Z})\mathrm{X})^{-1} \\
-(\mathrm{X}^{\mathrm{T}}(\mathrm{I} - \mathrm{M_Z})\mathrm{X})^{-1} & (\mathrm{X}^{\mathrm{T}}\mathrm{M_Z}\mathrm{X})^{-1} - (\mathrm{X}^{\mathrm{T}}(\mathrm{I} - \mathrm{M_Z})\mathrm{X})^{-1}
\end{bmatrix}
$$
代入上式知
$$
\begin{aligned}
\hat\beta_{\text{CF1}} &= ((\mathrm{X}^{\mathrm{T}}(\mathrm{I} - \mathrm{M_Z})\mathrm{X})^{-1}\mathrm{X}^{\mathrm{T}} - (\mathrm{X}^{\mathrm{T}}(\mathrm{I} - \mathrm{M_Z})\mathrm{X})^{-1}\mathrm{X}^{\mathrm{T}}\mathrm{M_Z})\mathrm{y} \\
&= (\mathrm{X}^{\mathrm{T}}(\mathrm{I} - \mathrm{M_Z})\mathrm{X})^{-1}\mathrm{X}^{\mathrm{T}}(\mathrm{I} - \mathrm{M_Z})\mathrm{y} \\
&= (\mathrm{X}^{\mathrm{T}}(\mathrm{I} - \mathrm{M_Z})(\mathrm{I} - \mathrm{M_Z})\mathrm{X})^{-1}\mathrm{X}^{\mathrm{T}}(\mathrm{I} - \mathrm{M_Z})\mathrm{y} \\
&= (\underbrace{\mathrm{X}^{\mathrm{T}}\mathrm{Z}(\mathrm{Z}^{\mathrm{T}}\mathrm{Z})^{-1}\mathrm{Z}^{\mathrm{T}}}
\underbrace{\mathrm{Z}(\mathrm{Z}^{\mathrm{T}}\mathrm{Z})^{-1}\mathrm{Z}^{\mathrm{T}}\mathrm{X}})^{-1}\underbrace{\mathrm{X}^{\mathrm{T}}\mathrm{Z}(\mathrm{Z}^{\mathrm{T}}\mathrm{Z})^{-1}\mathrm{Z}^{\mathrm{T}}}\mathrm{y} \\
&= (\hat{\mathrm{X}}^{\mathrm{T}}\hat{\mathrm{X}})^{-1}\hat{\mathrm{X}}^{\mathrm{T}}\mathrm{y} \\
&= \hat\beta_{\text{2SLS}}
\end{aligned}
$$
原命题得证，这里同时用到了矩阵$\mathrm{I} - \mathrm{M_Z}$的对称幂等性。

另外可以关注$\text{Var}(\hat\beta_{\text{CF}})$与$\text{Var}(\hat\beta_{\text{2SLS}})$。我们知道
$$
\begin{aligned}
\text{Var}(\hat\beta_{\text{CF}}) &= \sigma_{\text{CF}}^2 \begin{bmatrix}
\mathrm{X}^{\mathrm{T}}\mathrm{X} & \mathrm{X}^{\mathrm{T}}\mathrm{M_Z}\mathrm{X} \\
\mathrm{X}^{\mathrm{T}}\mathrm{M_Z}\mathrm{X} & \mathrm{X}^{\mathrm{T}}\mathrm{M_Z}\mathrm{X}
\end{bmatrix}^{-1} \\
&= \sigma_{\text{CF}}^2 \begin{bmatrix}
(\mathrm{X}^{\mathrm{T}}(\mathrm{I} - \mathrm{M_Z})\mathrm{X})^{-1} & -(\mathrm{X}^{\mathrm{T}}(\mathrm{I} - \mathrm{M_Z})\mathrm{X})^{-1} \\
-(\mathrm{X}^{\mathrm{T}}(\mathrm{I} - \mathrm{M_Z})\mathrm{X})^{-1} & (\mathrm{X}^{\mathrm{T}}\mathrm{M_Z}\mathrm{X})^{-1} - (\mathrm{X}^{\mathrm{T}}(\mathrm{I} - \mathrm{M_Z})\mathrm{X})^{-1}
\end{bmatrix}
\end{aligned}
$$
即$\text{Var}(\hat\beta_{\text{CF1}}) = \sigma_{\text{CF}}^2 (\mathrm{X}^{\mathrm{T}}(\mathrm{I} - \mathrm{M_Z})\mathrm{X})^{-1} = \sigma_{\text{CF}}^2 (\hat{\mathrm{X}}^{\mathrm{T}}\hat{\mathrm{X}})^{-1}$ 。而$\text{Var}(\hat\beta_{\text{2SLS}}) = \sigma_{\text{2SLS}}^2 (\hat{\mathrm{X}}^{\mathrm{T}}\hat{\mathrm{X}})^{-1}$，两者只相差一个常数。这说明，无论用$\text{CF}$还是$\text{2SLS}$，其变量对应的回归系数估计值都相同，同时任意量系数的协方差估计之比也相同，即：对$\forall 1 \le i, j, m, n \le K$，有
$$
\hat\beta_{\text{CF1},i} = \hat\beta_{\text{2SLS},i} \\
\frac{cov(\hat\beta_{\text{CF1},i}, \hat\beta_{\text{CF1},j})}{cov(\hat\beta_{\text{2SLS},i}, \hat\beta_{\text{2SLS},j})} = \frac{cov(\hat\beta_{\text{CF1},m}, \hat\beta_{\text{CF1},n})}{cov(\hat\beta_{\text{2SLS},m}, \hat\beta_{\text{2SLS},n})}
$$
