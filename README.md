# baby-kyber-rs

## Kyber
Kyber is a post-quantum public-key encryption system. 
Its main use case is to establish keys of symmetric-key systems in higher-level protocols like TLS, Signal or OpenPGP. it is a post-quantum system because Kyber is specifically designed to be secure even in the presence of quantum computers.

In 2021 NIST decided it is worthy of standardization.


### BabyKyber
BabyKyber is a down-scaled version of the system. 
It's equivalent to “regular” Kyber, except that the security parameters are smaller, making everything much more readable. 

Also, “regular” Kyber uses compression for ciphertexts, which we will omit here as it is not relevant for the underlying crypto system.



## MLWE
The trick is that it is a hard problem to recover `s` from `(A t)`. 
In fact, recovering `s` would require an attacker to solve the module-learning-with-errors (MLWE) problem, on which this system is built. 
The MLWE problem is expected to be hard even for quantum computers, which is why it is used in PQC.





# Reference
* https://cryptopedia.dev/posts/kyber/
