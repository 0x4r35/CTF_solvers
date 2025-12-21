#challenge 3 The Strongest Exponent I Thought
from Crypto.Util.number import long_to_bytes, inverse

n = 11577571088650843455101259498919170514888559114754991176367426310626451800000708972242631218783418858053014238827368727619121931192815809305245124163444094802225405870297050523228593276632591856050595253322599364103747724016666850166209913389873502963526134123874191503787054574544312952756389035649345032890497753865563100691095010626694207361877828026699086860070155460530090806390187874660145664212939747089918283596838097941086483857240175985051272705286675007014016041667499050742506806075203343665660277750845445942377591526363631000954388313070148480326021369420859362635929578502380138296606264683203983286197
e_val = 7295003680573243062191081850447296742887316782015364650489618985694064274472461768216735544329032076830044142827256877780040853871726633042483297304002890283242503674526130070158328305692362656491270378019397034874005287943151447421098854586070873187809805226478607730681630972237698922917494163140956836788
c = 4076499216416562896607441161849625990047067800665233407329196341540304308648851286217843088574720496112040934416510938682846567006808806010955921894978744639999560939385998309826612716572297961147343647566989002739744629380117872478800844703421167539558650736951573012109179090419063845190127599968982469602914208703636984004199567214853971463422671349587993404007992477629764085581293548006892140104708828058826600334161006488685079586296737011303994925944419846262719619542090760222527573288314177260270227503471544970415885712109366649795000097050560179546478161308265772556681319650087987273824051482087965825610

def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

def find_factors():
    candidates = [(0, 0)]
    for i in range(1024):
        next_candidates = []
        e_bit = (e_val >> i) & 1
        n_bit = (n >> i) & 1
        mask = (1 << (i + 1)) - 1

        for p_curr, q_curr in candidates:
            for p_bit in (0, 1):
                for q_bit in (0, 1):
                    if (p_bit ^ q_bit) != e_bit:
                        continue
                    p_new = p_curr | (p_bit << i)
                    q_new = q_curr | (q_bit << i)
                    if ((p_new * q_new) & mask) == (n & mask):
                        next_candidates.append((p_new, q_new))
        candidates = next_candidates
        if not candidates:
            return None

    for p_sol, q_sol in candidates:
        if p_sol * q_sol == n:
            return p_sol, q_sol
    return None

def solve_recursive(c_curr, e_curr, p_mod):
    g = gcd(e_curr, p_mod - 1)

    if g == 1:
        d = inverse(e_curr, p_mod - 1)
        return [pow(c_curr, d, p_mod)]

    if e_curr % 2 == 0:
        if pow(c_curr, (p_mod - 1) // 2, p_mod) != 1:
            return []

        sqrt_c = pow(c_curr, (p_mod + 1) // 4, p_mod)

        roots = []
        roots += solve_recursive(sqrt_c, e_curr // 2, p_mod)
        roots += solve_recursive(p_mod - sqrt_c, e_curr // 2, p_mod)
        return roots

    return []

def main():
    factors = find_factors()
    if not factors:
        return

    p, q = factors

    roots_p = solve_recursive(c, e_val, p)
    roots_q = solve_recursive(c, e_val, q)

    inv_p = inverse(p, q)
    inv_q = inverse(q, p)

    for rp in roots_p:
        for rq in roots_q:
            m = (rp * q * inv_q + rq * p * inv_p) % n
            try:
                flag_bytes = long_to_bytes(m)
                if b'TSGCTF' in flag_bytes:
                    print(flag_bytes.decode())
                    return
            except:
                continue

if __name__ == "__main__":
    main()
 
# FLAG :  TSGCTF{Wait,_caret_is_XOR,_not_power!}
