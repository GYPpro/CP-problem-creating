---
title: "我的题目"
author: "出题人"
date: "2025-04-23"
mainfont: "SimSun"            # 正文字体：系统里装了的字体名
sansfont: "Microsoft YaHei"   # 无衬线字体
monofont: "Fira Code"         # 等宽字体
geometry: "margin=2cm"        # （可选）页边距
---

# Problem. Embrace the Solitude of Destiny

$\begin{array}{l l}Time \; limit \; per \; test & :  1 \;  seconds\\ Memory \; limit \; per \; test &:  512 \;  megabytes\\Input& : standard\;  input \\Output & :standard\;  output \\ \end{array}$

> 塞纳托斯是死亡和湮灭的泰坦，无人知晓其真正的形象——因为死者无法将传说告知生者。据说。它的触碰会令花朵枯萎，土壤开裂，金属锈蚀，血肉成为灰烬——所有的一切被包裹在浓黑的雾气里，解离为没有意义的碎片。
>
> 黄金裔虽然自视为半神的血脉，但也无法逃脱有机生命应有的结局，死去。冥河会寻找死亡，并在空无一物的旷野上汇聚成通往冥界的道路。当战争夺去过多的生命，众多灵魂无处安息时，冥河的潮汐便会涌上大地，而塞纳托斯则会掌船前来，接引亡魂去往来世。人们认为，它的公正仅次于塔兰顿，因为它只取走众神手中多余之物。
> 
> 随着世界入夜，死亡与黑夜成为了翁法罗斯的主宰，未知的黑潮已卷席着巨大的浪涛吞没凡界，将冥界与现世的平衡打破。而在诸多战场上游荡的亡灵和怪物，无法流转的生死似乎已昭示，即使是掌管死亡的塞纳托斯，也对其无能为力。

作为哀地里亚的督战圣女，重新汇聚通往冥界的道路、寻索塞纳托斯的火种，是 Castorice 的责任。

![alt text](image.png)

$$
Figure : 「死荫的侍女」Castorice
$$


由于 Castorice 身负名为「死亡之触」的诅咒，“指尖触碰到的生命都会凋亡归无”，无论黄金裔或是元老院都无法靠近她，更遑论在向泰坦、向神格的求索过程中帮助她。而这，也一定程度上成为了 Castorice 命途上的阻碍。

事情的转机源自于 Stelle 的到来。这位来自翁法罗斯之外的少女不受「死亡之触」的影响，顺理成章地受黄金裔所托来与 Castorice 共同面对命运。

经过一段时间的交流与探索， Stelle 敏锐地发现，翁法罗斯的所谓「冥河所在的旷野」，本质是一个由向量化表述的，泰坦、半神、命途与四维时空张成的七维希尔伯特空间。那么， Castorice 的求索之路可以用一段游走来表示：

在希尔伯特空间 $\mathcal{H} = \mathcal{H}_C \otimes \mathcal{H}_P \otimes \mathcal{H}_R $ 上，规定条件置换 $U = S \cdot (C \otimes R)$ ，系统初态 $ |\psi_0\rangle = (a|L\rangle + b|R\rangle + c|U\rangle ) \otimes |\phi_{x_0}\rangle$，幺正算符$C$，那么作用在基矢上的置换算子可以表示为
$$
S|L,\phi_{x_0}\rangle = |\phi_{x_0 - r}\rangle,S|R,\phi_{x_0}\rangle = |L,\phi_{x_0 - r}\rangle,S|U,\phi_{x_0}\rangle = |U,\phi_{x_0 - r}\rangle
$$

显然，这些变换均为线性变换。

同时，由于非正交量子态不可区分，求索之路的成功有不少对节点的限制。那些条件的求解过于困难， Stelle 将其简化到了最基本的二维情况（只保留命途与时间两个维度），希望在求解的过程中找到更加一般的数学规律。

不过 Stelle 和 Castorice 都不是算法领域大神，难以应对在空间中规模庞大的限制条件数据。她们将问题上传到了JNU（Jupyter Neural-network Union，神经网络交互计算协会）进行悬赏，以求找到高效的解决方案。

在该问题中，Stelle 将 Castorice 的求索之路表述成在希尔伯特空间中的一连串置换、限制条件为一系列同余方程组，初始值则由 Castorice 和她自己组成的二维向量表示。特别地， Castorice 所代表的值已经被翁法罗斯的命运固定，这条路的一切变数都取决于天外来客 Stelle。

当然，在限制条件下，方程可能没有解，此时 Castorice 的故事就不可避免地走向悲剧的结局。

作为 JNU 的一员，你虽然不缺那点悬赏，但对这个问题本身非常感兴趣，因此决定稍作尝试。


**具体地：**

对于初始值向量 $|A_0\rangle = [x_0^{(1)},x_0^{(2)}]^T$ $^\ddagger$， Stelle 会给出 $K$ 个阶段的置换。其中，第 $i$ 个阶段的每次置换 $|A_{i,k-1}\rangle \mapsto |A_{i,k}\rangle$可以表示为：
$
\begin{cases}
x_{i,k}^{(1)} = a_i^{(1)} x_{i,k-1}^{(1)} + b_i^{(1)} x_{i,k-1}^{(2)} \\
x_{i,k}^{(2)} = a_i^{(2)} x_{i,k-1}^{(1)} + b_i^{(2)} x_{i,k-1}^{(2)}
\end{cases}
$

并重复执行 $r_i$ 次。

在该置换**结束**后，有当前值 $|A_{i,r_i}\rangle^T = [x_{i}^{(1)},x_{i}^{(2)}]$ ，需要满足：

$$
\begin{cases}
x_{i}^{(1)} \equiv c_i^{(1)} (mod \; m_i) \\
x_{i}^{(2)} \equiv c_i^{(2)} (mod \; m_i) 
\end{cases}
$$

现在，对于给定的 $K$ 个阶段的置换和其中一个初始值 $x_0^{(1)}$ ，你需要确定最小的非负整数初始值 $x_0^{(2)}$，使得所有给定限制条件成立。

$\ddagger$：此处的上角标为编号，并非乘方，下同。

<!-- **具体地，在本问题中**：将 Castorice 的求索之路表述成在欧几里得空间中一维初始值上的一串线性变换。 Stelle 会给出 $K$ 个阶段的线性变换，每个阶段的线性变换可以表示为 $  a_i x + b_i \mapsto x $ ，并重复执行 $r_i$ 次。

在第 $j$ 阶段结束后，当前的 $x$ 值需要满足 $x \equiv c_j ( mod \; m_j )$，才认为本段路径合法。

你需要确定一个可行的最小正整数 $x_0$ ，使得所有限制条件成立。

形式化地，有同余方程组 
$$
\begin{cases} 
x_{r_1} & \equiv c_1 ( mod \; m_1 ) \\
x_{r_1+r_2} & \equiv c_2 ( mod \; m_2 ) \\
& \vdots \\
x_{\sum_{i=1}^K r_i} & \equiv c_K ( mod \; m_K ) \\
\end{cases}
$$

其中 $\forall i \lt 1 , x_{i} = a_k x_{i-1} + b_k $ , $k$ 满足 $ \sum _{j=1}^{k-1} r_j < i \le \sum _{j=1}^{k} r_j $. -->

## Input

第一行两个整数 $K,x_0^{(1)}$， $1\le K \le 5\times 10^2,1\le x_0^{(1)} \le 1\times 10^9$，代表变换的段数和其中一个初始值。

后续，对于每段变换：包含三行输入。

第一行为一个整数 $r_i$ ，代表变换段数。

第二行为三个整数 $c_i^{(1)}, c_i^{(2)}, m_i$ ，含义如题面所示，代表置换后限制的同余方程参数。

第三行为四个整数 $ a_i^{(1)}, b_i^{(1)},  a_i^{(2)}, b_i^{(2)} $ ，含义如题面所示，代表置换参数。


其中，对于所有的 $a,b,c,d,r$，有 $1 \le a,b,c,d,r \le 1 \times 10^9$。

数据保证：所有 $m_i$ 的 LCM 不超过 $1\times 10^{18}$。

## Output

一行一个整数，代表 $x_0$。

不存在非负整数解时，输出 -1 。

### Input Example 1

```text
1 2
2
3 0 5
1 1 4 5
```

### Output Example 1

```text
3
```

<!-- 



### OUTGT

我们知道，一般的线性同余方程组，形如$x \equiv b_i (mod \; m_i)$，可以用中国剩余定理合并解系来求得最非负整数解。

对于另一形式的同余方程组，形如 $a_i \times x \equiv b_i (mod \; m_i)$，在$a_i$ 与 $m_i$ 不互质时不能通过逆元简化到一般的形式。那么是否有类似的合并解系的解法？ -->

# Problem. Embrace the Solitude of Destiny

$\begin{array}{l l}Time \; limit \; per \; test & :  1 \;  seconds\\ Memory \; limit \; per \; test &:  512 \;  megabytes\\Input& : standard\;  input \\Output & :standard\;  output \\ \end{array}$

> 塞纳托斯是死亡和湮灭的泰坦，无人知晓其真正的形象——因为死者无法将传说告知生者。据说。它的触碰会令花朵枯萎，土壤开裂，金属锈蚀，血肉成为灰烬——所有的一切被包裹在浓黑的雾气里，解离为没有意义的碎片。
>
> 黄金裔虽然自视为半神的血脉，但也无法逃脱有机生命应有的结局，死去。冥河会寻找死亡，并在空无一物的旷野上汇聚成通往冥界的道路。当战争夺去过多的生命，众多灵魂无处安息时，冥河的潮汐便会涌上大地，而塞纳托斯则会掌船前来，接引亡魂去往来世。人们认为，它的公正仅次于塔兰顿，因为它只取走众神手中多余之物。
> 
> 随着世界入夜，死亡与黑夜成为了翁法罗斯的主宰，未知的黑潮已卷席着巨大的浪涛吞没凡界，将冥界与现世的平衡打破。而在诸多战场上游荡的亡灵和怪物，无法流转的生死似乎已昭示，即使是掌管死亡的塞纳托斯，也对其无能为力。

作为哀地里亚的督战圣女，重新汇聚通往冥界的道路、寻索塞纳托斯的火种，是 Castorice 的责任。

![alt text](image.png)

$$
Figure : 「死荫的侍女」Castorice
$$


由于 Castorice 身负名为「死亡之触」的诅咒，“指尖触碰到的生命都会凋亡归无”，无论黄金裔或是元老院都无法靠近她，更遑论在向泰坦、向神格的求索过程中帮助她。而这，也一定程度上成为了 Castorice 命途上的阻碍。

事情的转机源自于 Stelle 的到来。这位来自翁法罗斯之外的少女不受「死亡之触」的影响，顺理成章地受黄金裔所托来与 Castorice 共同面对命运。

经过一段时间的交流与探索， Stelle 敏锐地发现，翁法罗斯的所谓「冥河所在的旷野」，本质是一个由向量化表述的，泰坦、半神、命途与四维时空张成的七维希尔伯特空间。那么， Castorice 的求索之路可以用一段游走来表示：

在希尔伯特空间 $\mathcal{H} = \mathcal{H}_C \otimes \mathcal{H}_P \otimes \mathcal{H}_R $ 上，规定条件置换 $U = S \cdot (C \otimes R)$ ，系统初态 $ |\psi_0\rangle = (a|L\rangle + b|R\rangle + c|U\rangle ) \otimes |\phi_{x_0}\rangle$，幺正算符$C$，那么作用在基矢上的置换算子可以表示为
$$
S|L,\phi_{x_0}\rangle = |\phi_{x_0 - r}\rangle,S|R,\phi_{x_0}\rangle = |L,\phi_{x_0 - r}\rangle,S|U,\phi_{x_0}\rangle = |U,\phi_{x_0 - r}\rangle
$$

显然，这些变换均为线性变换。

同时，由于非正交量子态不可区分，求索之路的成功有不少对节点的限制。那些条件的求解过于困难， Stelle 将其简化到了最基本的二维情况（只保留命途与时间两个维度），希望在求解的过程中找到更加一般的数学规律。

不过 Stelle 和 Castorice 都不是算法领域大神，难以应对在空间中规模庞大的限制条件数据。她们将问题上传到了JNU（Jupyter Neural-network Union，神经网络交互计算协会）进行悬赏，以求找到高效的解决方案。

在该问题中，Stelle 将 Castorice 的求索之路表述成在希尔伯特空间中的一连串置换、限制条件为一系列同余方程组，初始值则由 Castorice 和她自己组成的二维向量表示。特别地， Castorice 所代表的值已经被翁法罗斯的命运固定，这条路的一切变数都取决于天外来客 Stelle。

当然，在限制条件下，方程可能没有解，此时 Castorice 的故事就不可避免地走向悲剧的结局。

作为 JNU 的一员，你虽然不缺那点悬赏，但对这个问题本身非常感兴趣，因此决定稍作尝试。


**具体地：**

对于初始值向量 $|A_0\rangle = [x_0^{(1)},x_0^{(2)}]^T$ $^\ddagger$， Stelle 会给出 $K$ 个阶段的置换。其中，第 $i$ 个阶段的每次置换 $|A_{i,k-1}\rangle \mapsto |A_{i,k}\rangle$可以表示为：
$
\begin{cases}
x_{i,k}^{(1)} = a_i^{(1)} x_{i,k-1}^{(1)} + b_i^{(1)} x_{i,k-1}^{(2)} \\
x_{i,k}^{(2)} = a_i^{(2)} x_{i,k-1}^{(1)} + b_i^{(2)} x_{i,k-1}^{(2)}
\end{cases}
$

并重复执行 $r_i$ 次。

在该置换**结束**后，有当前值 $|A_{i,r_i}\rangle^T = [x_{i}^{(1)},x_{i}^{(2)}]$ ，需要满足：

$$
\begin{cases}
x_{i}^{(1)} \equiv c_i^{(1)} (mod \; m_i) \\
x_{i}^{(2)} \equiv c_i^{(2)} (mod \; m_i) 
\end{cases}
$$

现在，对于给定的 $K$ 个阶段的置换和其中一个初始值 $x_0^{(1)}$ ，你需要确定最小的非负整数初始值 $x_0^{(2)}$，使得所有给定限制条件成立。

$\ddagger$：此处的上角标为编号，并非乘方，下同。

<!-- **具体地，在本问题中**：将 Castorice 的求索之路表述成在欧几里得空间中一维初始值上的一串线性变换。 Stelle 会给出 $K$ 个阶段的线性变换，每个阶段的线性变换可以表示为 $  a_i x + b_i \mapsto x $ ，并重复执行 $r_i$ 次。

在第 $j$ 阶段结束后，当前的 $x$ 值需要满足 $x \equiv c_j ( mod \; m_j )$，才认为本段路径合法。

你需要确定一个可行的最小正整数 $x_0$ ，使得所有限制条件成立。

形式化地，有同余方程组 
$$
\begin{cases} 
x_{r_1} & \equiv c_1 ( mod \; m_1 ) \\
x_{r_1+r_2} & \equiv c_2 ( mod \; m_2 ) \\
& \vdots \\
x_{\sum_{i=1}^K r_i} & \equiv c_K ( mod \; m_K ) \\
\end{cases}
$$

其中 $\forall i \lt 1 , x_{i} = a_k x_{i-1} + b_k $ , $k$ 满足 $ \sum _{j=1}^{k-1} r_j < i \le \sum _{j=1}^{k} r_j $. -->

## Input

第一行两个整数 $K,x_0^{(1)}$， $1\le K \le 5\times 10^2,1\le x_0^{(1)} \le 1\times 10^9$，代表变换的段数和其中一个初始值。

后续，对于每段变换：包含三行输入。

第一行为一个整数 $r_i$ ，代表变换段数。

第二行为三个整数 $c_i^{(1)}, c_i^{(2)}, m_i$ ，含义如题面所示，代表置换后限制的同余方程参数。

第三行为四个整数 $ a_i^{(1)}, b_i^{(1)},  a_i^{(2)}, b_i^{(2)} $ ，含义如题面所示，代表置换参数。


其中，对于所有的 $a,b,c,d,r$，有 $1 \le a,b,c,d,r \le 1 \times 10^9$。

数据保证：所有 $m_i$ 的 LCM 不超过 $1\times 10^{18}$。

## Output

一行一个整数，代表 $x_0$。

不存在非负整数解时，输出 -1 。

### Input Example 1

```text
1 2
2
3 0 5
1 1 4 5
```

### Output Example 1

```text
3
```

<!-- 



### OUTGT

我们知道，一般的线性同余方程组，形如$x \equiv b_i (mod \; m_i)$，可以用中国剩余定理合并解系来求得最非负整数解。

对于另一形式的同余方程组，形如 $a_i \times x \equiv b_i (mod \; m_i)$，在$a_i$ 与 $m_i$ 不互质时不能通过逆元简化到一般的形式。那么是否有类似的合并解系的解法？ -->