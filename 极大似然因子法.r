# S是样本方差阵或相关矩阵
# m是主因子的个数
# d是特殊方差的估计值
# output

## 极大似然因子分析
mle.fa=function(S, m, d){
    p=nrow(S)
    diag_S=diag(S)
    sum_rank=sum(diag_S)
    rowname=paste("X", 1:p, sep="")
    colname=paste("Factor", 1:m, sep="")
    A=matrix(0, nrow=p, ncol=m, 
            dimnames=list(rowname, colname))
    
    # kmax最大迭代次数
    kmax=20
    k=1
    repeat{
    d1=d; d2=1/sqrt(d); 
    eig=eigen(S * (d2 %o% d2))
    for (i in 1:m)
        A[,i]=sqrt(eig$values[i]-1)*eig$vectors[,i]
    A=diag(sqrt(d)) %*% A
    d=diag(S-A%*%t(A))
    if ((sqrt(sum((d-d1)^2))<1e-5)|k==kmax) break
    k=k+1
  }
  
  rowname=c("SS loadings","Proportion Var","Cumulative Var")
  f=matrix(0, nrow=3, ncol=m, dimnames=list(rowname, colname))
  for (i in 1:m){
    f[1,i]=sum(A[,i]^2)
    f[2,i]=f[1,i]/sum_rank
    f[3,i]=sum(f[1,1:i])/sum_rank
  }
  method=c("MLE")
  list(method=method, loadings=A, 
       var=cbind(common=diag_S-d, spcific=d),f=f,iterative=k) 
}

# 导入数据
setwd("~/Multi-Statis/hw4")
data = read.csv("Algerian_forest_fires_dataset_UPDATE.csv")
sdat = data[, 4:13]
sdat = as.data.frame(lapply(sdat,as.factor))
sdat = as.data.frame(lapply(sdat,as.numeric))
# 数据标准化
data = scale(sdat, scale=TRUE)
# data = iris[,1:4]
R = cor(data)

## 特殊方差d的选取
d = 1 / diag(solve(R)) 

## 选取因子个数为4
## 输出loadings-因子载荷，var-共性方差和特殊方差
## B-因子F对变量X的贡献：贡献率和累积贡献率，iterative-迭代次数
fa = mle.fa(R, m=4, d)
fa

# 近似S的误差平方和Q(m)
E = R-fa$loadings %*% t(fa$loadings)-diag(fa$var[,2])
print(sum(E^2))

# 检验因子之间的相关性
chisq.test(fa$f)  # p值为0.96 不相关


