from collections import Counter

AMINO_ACID_MASSES = [57,71,87,97,99,101,103,113,114,115,128,129,131,137,147,156,163,186]

def linear_spectrum(peptide):
    prefix=[0]
    for m in peptide:
        prefix.append(prefix[-1]+m)
    spec=[0]
    n=len(peptide)
    for i in range(n):
        for j in range(i+1,n+1):
            spec.append(prefix[j]-prefix[i])
    return sorted(spec)

def cyclic_spectrum(peptide):
    prefix=[0]
    for m in peptide: prefix.append(prefix[-1]+m)
    total=prefix[-1]
    spec=[0]
    n=len(peptide)
    for i in range(n):
        for j in range(i+1,n+1):
            sub=prefix[j]-prefix[i]
            spec.append(sub)
            if i>0 and j<n:
                spec.append(total-sub)
    return sorted(spec)

def score_spectrum(theo, exp_counter):
    ct=Counter(theo)
    return sum(min(ct[m], exp_counter.get(m,0)) for m in ct)

def expand(peptides):
    res=[]
    for pep,m in peptides:
        for a in AMINO_ACID_MASSES:
            res.append((pep+(a,), m+a))
    return res

def trim(peptides, exp_counter, N):
    scored=[(score_spectrum(linear_spectrum(list(pep)), exp_counter), pep, mass) for pep,mass in peptides]
    scored.sort(reverse=True,key=lambda x:x[0])
    if len(scored)<=N:
        return [(p,m) for s,p,m in scored]
    cutoff=scored[N-1][0]
    return [(p,m) for s,p,m in scored if s>=cutoff]

def leaderboard_all(spectrum, N=1000):
    exp_sorted=sorted(spectrum); exp_counter=Counter(exp_sorted)
    parent_mass=max(spectrum)
    leaderboard=[((),0)]
    best_score=-1; leaders=[]
    while leaderboard:
        leaderboard=expand(leaderboard)
        # dedup
        seen=set(); new=[]
        for pep,m in leaderboard:
            if pep in seen: continue
            seen.add(pep); new.append((pep,m))
        leaderboard=new
        new_board=[]
        for pep,m in leaderboard:
            if m==parent_mass:
                sc=score_spectrum(cyclic_spectrum(list(pep)), exp_counter)
                if sc>best_score:
                    best_score=sc; leaders=[pep]
                elif sc==best_score:
                    leaders.append(pep)
                new_board.append((pep,m))
            elif m<parent_mass:
                new_board.append((pep,m))
        leaderboard=trim(new_board, exp_counter, N)
    return leaders, best_score

# put your Spectrum25 here (or read from stdin/file)
spectrum = list(map(int, "0 97 99 113 114 115 128 128 147 147 163 186 227 241 242 244 244 256 260 261 262 283 291 309 330 333 340 347 385 388 389 390 390 405 435 447 485 487 503 504 518 544 552 575 577 584 599 608 631 632 650 651 653 672 690 691 717 738 745 770 779 804 818 819 827 835 837 875 892 892 917 932 932 933 934 965 982 989 1039 1060 1062 1078 1080 1081 1095 1136 1159 1175 1175 1194 1194 1208 1209 1223 1322".split()))
leaders, best = leaderboard_all(spectrum, N=1000)
formatted = ["-".join(map(str,p)) for p in leaders]
print("\n".join(formatted))
print(best)