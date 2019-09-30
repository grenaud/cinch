import numpy as np
from sklearn.decomposition import NMF, non_negative_factorization

R = [
     [5,3,0,1],
     [4,0,0,1],
     [1,1,0,5],
     [1,0,0,4],
     [0,1,5,4],
    ]

R = np.array(R)
nmf = NMF(n_components=2, init='random', random_state=0)

W = nmf.fit_transform(R);
H = nmf.components_;
nR = np.dot(W,H)

print(nmf)
print("W");
print(W);
print("H");
print(H);
print(R);
print(nR);
S = [
     [6,4,1,2],
    ]

nS=nmf.transform(S);
print(S);
print(nS);



fixed_H=H;
W, H, n_iter = non_negative_factorization(R, n_components=2, init='custom', random_state=0, update_H=False, H=fixed_H)


print("W");
print(W);
print("H");
print(H);
