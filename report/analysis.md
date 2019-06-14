# DNA Concatente
## Thoughts

* 观察以下交叠区域的平均长度，也可以通过复制次数进行理论计算。
* 感觉和n-gram问题有点相像，而且这里只有四个元素，一级词表很小。$ 4_^9 < 10^6 $ 也就是说词表可以大概存序列长度小于9的所有序列，这个完全存的下，而且也不会有什么问题。如果再涨以下，到10^12，有点大，这样就有18的长度，好像其实也没长多少。
* 每一个数据集的基因都是由五个物种拆出来的，也就是说有5个物种，但是物种之间的很有可能会share同样的基因片段，这样不能暴力去重了，找路的时候，找到一条最长路之后，可能还会有经过某一个路的可能。

在data1中，原则上如果有分叉，那么一定是原基因序列出现了重复片段
那么必然有环出现，且这个环实际上要走两遍才能拼出原序列。


## Reference
* https://academic.oup.com/nar/article/45/6/e43/2638393?sid=4c8c5d1a-1fe3-4702-a19e-1848b4c261d5 

de novo
http://www.sohu.com/a/202971027_278730
问一下同一个数据集内的五个物种时什么关系，不同数据集之间又时什么关系。序列切割时选用了什么算法，随机性体现在哪里（方便进行理论分析）

## 数学分析

$$
\text{Assume that the length of the genome is N, } \\
\text{and let k reprensent the length of kmers,} \\
\text{Then the number of kmers } N_{kmers} = N-k+1 \\
\text{The all possible value that a k-length kmer can be is } R_{kmers} = 4^k \\
P = \frac{N_{kmers} }{R_{kmers}} = \frac{N-k+1}{4^k} \\
\lim_{k>>\log_4^N} P = 0  \\
\text{Let's set N to 100000} \\
P=\{ \begin{align} 
& \le 0.5 &,k>=9 \\
& ... &, k<9
\end{align} \\
\text{This shows that the kmers is not likely to occur many times,}\\
\text{and there will be many kmers won't occur after cutting.  }
$$



[1]: https://alexbowe.com/succinct-debruijn-graphs/#our-repr

