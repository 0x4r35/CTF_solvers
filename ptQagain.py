#challenge 2 PTQagain 
import math

n = "given"
e = "given"
c = "given"
# Attack
# e = 1 mod (p-1) => e-1 is multiple of p-1
# 2^(e-1) = 1 mod p
# p = gcd(2^(e-1) - 1, n)
val = pow(2, e - 1, n) - 1

p = math.gcd(val, n)

q = n // p

phi = (p - 1) * (q - 1)

d = pow(e, -1, phi)

m = pow(c, d, n)

print(m.to_bytes((m.bit_length() + 7) // 8, 'big').decode())

#FLAG: TSGCTF{My_power_is_too_dangerous_to_use!} 