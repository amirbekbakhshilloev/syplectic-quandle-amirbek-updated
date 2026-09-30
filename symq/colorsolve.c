#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

typedef uint64_t word;
static int N, n, m, W;
static int *op;
static int (*rel)[3];
static int **by_arc; static int *by_cnt;
static word fullw[64];
static long long nsol = 0;
static int *queue; static char *inq;
static int *tmpx;
static int print_solutions = 1;
static int nsel = 0; static int sel[64];
static int npri = 0; static int pri[64];

static inline int popc(const word *d) { int s = 0; for (int i = 0; i < W; i++) s += __builtin_popcountll(d[i]); return s; }
static inline int isfull(const word *d) { for (int i = 0; i < W; i++) if (d[i] != fullw[i]) return 0; return 1; }
static inline int isempty(const word *d) { for (int i = 0; i < W; i++) if (d[i]) return 0; return 1; }
static inline int has(const word *d, int z) { return (d[z >> 6] >> (z & 63)) & 1; }
static inline void setb(word *d, int z) { d[z >> 6] |= ((word)1) << (z & 63); }

static int revise(int r, word *dom) {
    int a = rel[r][0], b = rel[r][1], c = rel[r][2];
    word *Da = dom + (size_t)a * W, *Db = dom + (size_t)b * W, *Dc = dom + (size_t)c * W;
    int fulls = isfull(Da) + isfull(Db) + isfull(Dc);
    if (fulls >= 2) return 0;
    word nA[W], nB[W], nC[W];
    memset(nA, 0, sizeof(word) * W); memset(nB, 0, sizeof(word) * W); memset(nC, 0, sizeof(word) * W);
    int la = 0;
    for (int i = 0; i < W; i++) { word t = Da[i]; while (t) { int bit = __builtin_ctzll(t); tmpx[la++] = i * 64 + bit; t &= t - 1; } }
    for (int i = 0; i < W; i++) {
        word t = Db[i];
        while (t) {
            int y = i * 64 + __builtin_ctzll(t); t &= t - 1;
            int hit = 0;
            for (int j = 0; j < la; j++) {
                int x = tmpx[j]; int z = op[x * N + y];
                if (has(Dc, z)) { setb(nC, z); setb(nA, x); hit = 1; }
            }
            if (hit) setb(nB, y);
        }
    }
    int ch = 0;
    if (memcmp(nA, Da, sizeof(word) * W)) { memcpy(Da, nA, sizeof(word) * W); ch |= 1; if (isempty(Da)) return -1; }
    if (memcmp(nB, Db, sizeof(word) * W)) { memcpy(Db, nB, sizeof(word) * W); ch |= 2; if (isempty(Db)) return -1; }
    if (memcmp(nC, Dc, sizeof(word) * W)) { memcpy(Dc, nC, sizeof(word) * W); ch |= 4; if (isempty(Dc)) return -1; }
    return ch;
}

static int propagate(word *dom, int *init, int ninit) {
    int qh = 0, qt = 0;
    memset(inq, 0, m);
    for (int i = 0; i < ninit; i++) { int r = init[i]; if (!inq[r]) { inq[r] = 1; queue[qt++] = r; } }
    while (qh != qt) {
        int r = queue[qh]; qh = (qh + 1) % (m + 1); inq[r] = 0;
        int ch = revise(r, dom);
        if (ch < 0) return 0;
        for (int role = 0; role < 3; role++) if (ch & (1 << role)) {
            int x = rel[r][role];
            for (int k = 0; k < by_cnt[x]; k++) { int rj = by_arc[x][k]; if (!inq[rj]) { inq[rj] = 1; queue[qt] = rj; qt = (qt + 1) % (m + 1); } }
        }
    }
    return 1;
}

static void dfs(word *dom) {
    int best = -1, bs = N + 1;
    for (int t = 0; t < npri; t++) { int i = pri[t]; int s = popc(dom + (size_t)i * W); if (s > 1) { best = i; bs = s; break; } }
    if (best < 0) for (int i = 0; i < n; i++) { int s = popc(dom + (size_t)i * W); if (s > 1 && s < bs) { bs = s; best = i; if (s == 2) break; } }
    if (best < 0) {
        int sol[n];
        for (int i = 0; i < n; i++) { const word *d = dom + (size_t)i * W; for (int j = 0; j < W; j++) if (d[j]) { sol[i] = j * 64 + __builtin_ctzll(d[j]); break; } }
        for (int r = 0; r < m; r++) if (op[sol[rel[r][0]] * N + sol[rel[r][1]]] != sol[rel[r][2]]) return;
        nsol++;
        if (print_solutions) {
            if (nsel) { for (int i = 0; i < nsel; i++) printf(i ? " %d" : "%d", sol[sel[i]]); printf("\n"); }
            else { for (int i = 0; i < n; i++) printf(i ? " %d" : "%d", sol[i]); printf("\n"); }
        }
        return;
    }
    size_t sz = (size_t)n * W;
    word *d2 = malloc(sizeof(word) * sz);
    word *Db = dom + (size_t)best * W;
    for (int j = 0; j < W; j++) {
        word t = Db[j];
        while (t) {
            int v = j * 64 + __builtin_ctzll(t); t &= t - 1;
            memcpy(d2, dom, sizeof(word) * sz);
            memset(d2 + (size_t)best * W, 0, sizeof(word) * W); setb(d2 + (size_t)best * W, v);
            if (propagate(d2, by_arc[best], by_cnt[best])) dfs(d2);
        }
    }
    free(d2);
}

int main(int argc, char **argv) {
    for (int ai = 1; ai < argc; ai++) {
        if (!strcmp(argv[ai], "-c")) print_solutions = 0;
        else if (!strcmp(argv[ai], "-b")) { npri = atoi(argv[++ai]); for (int i = 0; i < npri; i++) pri[i] = atoi(argv[++ai]); }
        else if (!strcmp(argv[ai], "-s")) { nsel = atoi(argv[++ai]); for (int i = 0; i < nsel; i++) sel[i] = atoi(argv[++ai]); }
    }
    if (scanf("%d %d %d", &N, &n, &m) != 3) return 1;
    W = (N + 63) / 64;
    for (int i = 0; i < W; i++) { int lo = i * 64, hi = (i + 1) * 64 < N ? (i + 1) * 64 : N; int k = hi - lo; fullw[i] = k == 64 ? ~(word)0 : ((((word)1) << k) - 1); }
    op = malloc(sizeof(int) * N * N);
    for (int i = 0; i < N * N; i++) if (scanf("%d", &op[i]) != 1) return 1;
    rel = malloc(sizeof(int[3]) * m);
    by_cnt = calloc(n, sizeof(int)); by_arc = malloc(sizeof(int *) * n);
    for (int i = 0; i < n; i++) by_arc[i] = malloc(sizeof(int) * 3 * m);
    for (int r = 0; r < m; r++) {
        if (scanf("%d %d %d", &rel[r][0], &rel[r][1], &rel[r][2]) != 3) return 1;
        for (int role = 0; role < 3; role++) {
            int x = rel[r][role], dup = 0;
            for (int k = 0; k < by_cnt[x]; k++) if (by_arc[x][k] == r) dup = 1;
            if (!dup) by_arc[x][by_cnt[x]++] = r;
        }
    }
    word *dom = malloc(sizeof(word) * (size_t)n * W);
    for (int i = 0; i < n; i++) memcpy(dom + (size_t)i * W, fullw, sizeof(word) * W);
    int k; if (scanf("%d", &k) == 1) {
        for (int t = 0; t < k; t++) {
            int arc, cnt; if (scanf("%d %d", &arc, &cnt) != 2) return 1;
            word *d = dom + (size_t)arc * W; memset(d, 0, sizeof(word) * W);
            for (int j = 0; j < cnt; j++) { int v; if (scanf("%d", &v) != 1) return 1; setb(d, v); }
        }
    }
    queue = malloc(sizeof(int) * (m + 1)); inq = calloc(m, 1); tmpx = malloc(sizeof(int) * N);
    int *all = malloc(sizeof(int) * m); for (int r = 0; r < m; r++) all[r] = r;
    if (propagate(dom, all, m)) dfs(dom);
    printf("# %lld\n", nsol);
    return 0;
}
