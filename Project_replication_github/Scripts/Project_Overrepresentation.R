#Hypergeometric test Mitochondrial proteins

# X is the number of positive proteins in your network, so the number of mitochondrial (or nuclear) proteins

# M is the total number of positive proteins, so for the mitochondrial ones that is the total number of mitochondrial proteins according to mitocarta

# N is the total number of negative proteins, so all proteins that are not mitochondrially located. According to the mitocarta list there are 19244 human genes and 19244-1158 =18086 so this is the number for the mitochondrial ones, for the nuclear ones you would need to substract the number of nuclear genes from 19244
#19244 total genes
#5916 nuclear genes
#1158 mitochondrial genes

# K is the sample size, which for our analyses is always 100.

x = 147
m = 5916
n = 13328
k = 200
dhyper(x = x, m = m, n = n, k = k)

#Expected number of mito proteins
k * m / (m + n)
