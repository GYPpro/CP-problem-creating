$$
线段 (x_1,y_1) (x_2,y_2) ，点(x_0,y_0)
先化作一般式：\\
A = y_2 - y_1, B = x_1 - x_2, C = x_2 y_1 - x_1 y_2\\
\\
t = \frac{Ax_0 + By_0 + C}{A^2 + B^2}\\
\\
那么  x_3 = x_0 + t \cdot (A,B)  y_3 = y_0 + t \cdot (A,B)\\
\\
对称点为 (x_3 * 2 - x_0, y_3 * 2 - y_0)\\
$$

# 1

对 $ u(t_{n+1}) $ 和 $ u'(t_{n+1}) $ 展开到 $\Delta t^3$：

$
u(t_{n+1}) = u(t_n) + \Delta t u'(t_n) + \frac{\Delta t^2}{2} u''(t_n) + \frac{\Delta t^3}{6} u'''(t_n) + O(\Delta t^4)
$  

$
u'(t_{n+1}) = u'(t_n) + \Delta t u''(t_n) + \frac{\Delta t^2}{2} u'''(t_n) + O(\Delta t^3)
$  

则有截断误差：  
$
\begin{align}
T_n  & \\
& = \frac{u(t_{n+1}) - u(t_n)}{\Delta t} - \frac{1}{2} [ u'(t_n) + u'(t_{n+1}) ] \\
 & = -\frac{\Delta t^2}{12} u'''(t_n) + O(\Delta t^3)  
\end{align}
$ 

精度阶数为2.

# 2
$
\begin{cases}
x'' = -\dfrac{k x}{(x^2 + y^2)^3}, \\
y'' = -\dfrac{k y}{(x^2 + y^2)^3}.
\end{cases}  
$ 

代入 $ u_1 = x, \, u_2 = x', \, u_3 = y, \, u_4 = y' $ 得
$
\begin{cases}
\dfrac{du_1}{dt} = u_2, \\
\dfrac{du_2}{dt} = -\dfrac{k u_1}{(u_1^2 + u_3^2)^3}, \\
\dfrac{du_3}{dt} = u_4, \\
\dfrac{du_4}{dt} = -\dfrac{k u_3}{(u_1^2 + u_3^2)^3}.
\end{cases}  
$  
则有迭代公式

$
\begin{cases}
u_1^{n+1} = u_1^n + \Delta t \cdot u_2^n, \\
u_2^{n+1} = u_2^n + \Delta t \cdot \left( -\dfrac{k u_1^n}{(u_1^n)^2 + (u_3^n)^2)^3} \right), \\
u_3^{n+1} = u_3^n + \Delta t \cdot u_4^n, \\
u_4^{n+1} = u_4^n + \Delta t \cdot \left( -\dfrac{k u_3^n}{(u_1^n)^2 + (u_3^n)^2)^3} \right).
\end{cases}  
$  
