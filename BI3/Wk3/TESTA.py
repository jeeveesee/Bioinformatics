def multiple_lcs(x, y, z):
    lx, ly, lz = len(x), len(y), len(z)

    # Build 3D DP table
    # dp[i][j][k] = LCS score for x[:i], y[:j], z[:k]
    dp = [[[0] * (lz + 1) for _ in range(ly + 1)] for _ in range(lx + 1)]

    for i in range(1, lx + 1):
        for j in range(1, ly + 1):
            for k in range(1, lz + 1):
                # All 7 transitions
                candidates = [
                    dp[i-1][j][k],       # advance x only (gap in y, z)
                    dp[i][j-1][k],       # advance y only
                    dp[i][j][k-1],       # advance z only
                    dp[i-1][j-1][k],     # advance x, y
                    dp[i-1][j][k-1],     # advance x, z
                    dp[i][j-1][k-1],     # advance y, z
                    dp[i-1][j-1][k-1] + (1 if x[i-1] == y[j-1] == z[k-1] else 0),  # advance all
                ]
                dp[i][j][k] = max(candidates)

    score = dp[lx][ly][lz]

    # Backtrack to find alignment
    ax, ay, az = [], [], []
    i, j, k = lx, ly, lz

    while i > 0 or j > 0 or k > 0:
        # Determine which transition led to current cell
        cur = dp[i][j][k]

        if i > 0 and j > 0 and k > 0 and cur == dp[i-1][j-1][k-1] + (1 if x[i-1] == y[j-1] == z[k-1] else 0):
            ax.append(x[i-1]); ay.append(y[j-1]); az.append(z[k-1])
            i -= 1; j -= 1; k -= 1
        elif i > 0 and j > 0 and cur == dp[i-1][j-1][k]:
            ax.append(x[i-1]); ay.append(y[j-1]); az.append('-')
            i -= 1; j -= 1
        elif i > 0 and k > 0 and cur == dp[i-1][j][k-1]:
            ax.append(x[i-1]); ay.append('-'); az.append(z[k-1])
            i -= 1; k -= 1
        elif j > 0 and k > 0 and cur == dp[i][j-1][k-1]:
            ax.append('-'); ay.append(y[j-1]); az.append(z[k-1])
            j -= 1; k -= 1
        elif i > 0 and cur == dp[i-1][j][k]:
            ax.append(x[i-1]); ay.append('-'); az.append('-')
            i -= 1
        elif j > 0 and cur == dp[i][j-1][k]:
            ax.append('-'); ay.append(y[j-1]); az.append('-')
            j -= 1
        else:
            ax.append('-'); ay.append('-'); az.append(z[k-1])
            k -= 1

    # Reverse since we backtracked
    ax = ''.join(reversed(ax))
    ay = ''.join(reversed(ay))
    az = ''.join(reversed(az))

    return score, ax, ay, az


if __name__ == "__main__":
    x = "ATATCCG"
    y = "TCCGA"
    z = "ATGTACTG"

    score, ax, ay, az = multiple_lcs(x, y, z)
    print(f"Score: {score}")
    print(ax)
    print(ay)
    print(az)

    # Verify: count columns where all three are same non-gap character
    verified_score = sum(
        1 for a, b, c in zip(ax, ay, az)
        if a == b == c and a != '-'
    )
    print(f"\nVerified score from alignment: {verified_score}")
    print(f"Alignment lengths match: {len(ax) == len(ay) == len(az)}")

    # Verify original strings are preserved (ignoring gaps)
    print(f"x preserved: {''.join(c for c in ax if c != '-') == x}")
    print(f"y preserved: {''.join(c for c in ay if c != '-') == y}")
    print(f"z preserved: {''.join(c for c in az if c != '-') == z}")