#challenge 1 bingo
from pwn import *

context.log_level = 'debug'

def solve():
    r = remote('35.194.98.181', 10961)
    
    r.recvuntil(b'input message (hex): ')

    payload_hex = b"Get Flag.".hex()
    r.sendline(payload_hex.encode())

    r.recvuntil(b'input sinature (int):')
    
    r.sendline(b'1337')

    r.interactive()

    if __name__ == '__main__' :
        solve()

 #FLAG:   TSGCTF{+he|23_was_4_f4rm3|2_had_a_|)o6_a|^_||)_8i|^_|6o_wa5_his_nam3}     
    