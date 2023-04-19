#include <bits/stdc++.h>
#include <sys/time.h>
#include <atcoder/all>

using namespace std;

#pragma region prototype_declaration
/* ============================================== 
    プロトタイプ宣言はここから
   ============================================== */

/*乱数生成器*/
struct RandGenerator {
    random_device seed_gen;
    mt19937 engine;
    mt19937_64 engine64;
    static const int pshift = 1000000000;
    RandGenerator() : engine(seed_gen()), engine64(seed_gen()) {}
    /*mod以下の乱数を返す（32bit）*/
    int rand(int mod) {
        return engine() % mod;
    }
    /*mod以下の乱数を返す（64bit）*/
    long long randll(long long mod) {
        return engine64() % mod;
    } 
    /*確率pでTrueを返す*/
    bool pjudge(double p) {
        int p_int;
        if(p > 1) p_int = pshift;
        else p_int = p * pshift;
        return rand(pshift) < p_int;
    }
} ryuka;

/*タイマー*/
struct Timer {
    double global_start;
    /*現在の時刻を返す*/
    double gettime() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return tv.tv_sec + tv.tv_usec * 1e-6;
    }
    void init() {
        global_start = gettime();
    }
    /*プログラム開始からの経過時間を返す*/
    double elapsed() {
        return gettime() - global_start;
    }
} toki;

struct Input {
    /*TODO: ここに入力変数を定義する*/
    const int n = 1000;
    vector<int> a, b, c, d;
    void read();
    pair<int,int> src(int id) const;
    pair<int,int> dst(int id) const;
} input;

struct Node {
    pair<int, int> pos;
    int p, id;
    Node(pair<int,int> pos, int p, int id) : pos(pos), p(p), id(id) {}
};

struct Output {
    /*TODO: ここに出力変数を定義する*/
    vector<int> nodes;
    vector<Node> route;
    vector<vector<int>> route_id;
    map<pair<int,int>,int> node_counter;
    Output();
    void print();
};

/*解を管理するクラス*/
struct State {
    Output output;
    long long score;
    long long length;
    State() : score(0) {}
    static State initState();
    static State initState(const vector<int>& nodes);
    static State generateState(const State& input_state);
};

struct State2 : State {
    static State2 generateState(const State2& input_state);
    static State2 convert(const State& input_state);
};

/*イテレーション管理クラス*/
template<class STATE>
struct IterationControl {
    int iteration_counter;
    int swap_counter;
    double average_time;
    double start_time;
    IterationControl() : iteration_counter(0), swap_counter(0) {}
    /*山登り法*/
    STATE climb(double time_limit, STATE initial_state) {
        start_time = toki.gettime();
        average_time = 0;
        STATE best_state = initial_state;
        double time_stamp = start_time;
        cerr << "[INFO] - IterationControl::climb - Starts climbing...\n";
        while(time_stamp - start_time + average_time < time_limit) {
            STATE current_state = STATE::generateState(best_state);
            if(current_state.score > best_state.score) {
                swap(best_state, current_state);
                swap_counter++;
            }
            iteration_counter++;
            time_stamp = toki.gettime();
            average_time = (time_stamp - start_time) / iteration_counter;
        }
        cerr << "[INFO] - IterationControl::climb - Iterated " << iteration_counter << " times and swapped " << swap_counter << " times.\n";
        return best_state;
    }
    /*焼きなまし法*/
    STATE anneal(double time_limit, double temp_start, double temp_end, STATE initial_state) {
        start_time = toki.gettime();
        average_time = 0;
        STATE best_state = initial_state;
        double elapsed_time = 0;
        cerr << "[INFO] - IterationControl::anneal - Starts annealing...\n";
        while(elapsed_time + average_time < time_limit) {
            double normalized_time = elapsed_time / time_limit;
            double temp_current = pow(temp_start, 1.0 - normalized_time) * pow(temp_end, normalized_time);
            STATE current_state = STATE::generateState(best_state);
            long long delta = current_state.score - best_state.score;
            if(delta > 0 || ryuka.pjudge(exp(1.0 * delta / temp_current)) ) {
                swap(best_state, current_state);
                swap_counter++;
            }
            iteration_counter++;
            elapsed_time = toki.gettime() - start_time;
            average_time = elapsed_time / iteration_counter;
        }
        cerr << "[INFO] - IterationControl::anneal - Iterated " << iteration_counter << " times and swapped " << swap_counter << " times.\n";
        return best_state;
    }
};

namespace Utils {
    long long calcLength(const Output& output);
    long long calcScore(const Output& output);
    long long calcScore(long long length);
    int calcDist(const pair<int,int>& a, const pair<int,int>& b);
    Output runInsertTSP(const vector<int>& nodes);
};

/* ============================================== 
    プロトタイプ宣言はここまで
   ============================================== */

#pragma endregion prototype_declaration

/*TODO: ここで入力を受け取る*/
void Input::read() {
    a.resize(n);
    b.resize(n);
    c.resize(n);
    d.resize(n);
    for(int i = 0; i < n; i++) {
        cin >> a[i] >> b[i] >> c[i] >> d[i];
    }
}

pair<int,int> Input::src(int id) const {
    return {a[id], b[id]};
}

pair<int,int> Input::dst(int id) const {
    return {c[id], d[id]};
}

/*TODO：ここで出力変数を初期化する。vectorのメモリ確保など*/
Output::Output() {
}

/*TODO：ここで答えを出力する*/
void Output::print() {
    cout << nodes.size();
    for(int id: nodes) cout << " " << id + 1;
    cout << endl;
    cout << route.size();
    for(Node node: route) cout << " " << node.pos.first << " " << node.pos.second;
    cout << endl;
}

/*TODO: ここで初期解を作成する*/
State State::initState() {
    vector<int> nodes;
    // (400,400)に近い方から50個選ぶ
    vector<int> idx(input.n);
    iota(idx.begin(), idx.end(), 0);
    sort(idx.begin(), idx.end(), [&](int l, int r){
        int ldist = max(Utils::calcDist({400, 400}, input.src(l)), Utils::calcDist({400, 400}, input.dst(l)));
        int rdist = max(Utils::calcDist({400, 400}, input.src(r)), Utils::calcDist({400, 400}, input.dst(r)));
        return ldist < rdist; 
    });
    for(int i = 0; i < 50; i++) nodes.push_back(idx[i]);
    return initState(nodes);
}

State State::initState(const vector<int>& nodes) {
    // とりあえず、挿入法で適当な巡回経路を作る。
    State res;
    res.output = Utils::runInsertTSP(nodes);
    res.length = Utils::calcLength(res.output);
    res.score = Utils::calcScore(res.output);
    return res;
}

/*TODO: ここでinput_stateを変化させた解を作る（局所探索）*/
State2 State2::generateState(const State2& input_state) {
    State2 res = input_state;
    // 2-opt
    int a = ryuka.rand(res.output.route.size() - 2) + 1;
    int b = ryuka.rand(res.output.route.size() - 2) + 1;
    auto check = [&](int x, int y) -> bool {
        Node v = res.output.route[x];
        assert(v.id >= 0);
        if(v.p == 0) {
            if(res.output.route_id[1][v.id] <= y) return false; 
        } 
        if(v.p == 1) {
            if(res.output.route_id[0][v.id] >= y) return false;
        }
        return true;
    };
    if(check(a, b) && check(b, a)) {
        Node va = res.output.route[a];
        Node vb = res.output.route[b];
        res.output.route_id[va.p][va.id] = b;
        res.output.route_id[vb.p][vb.id] = a;
        res.length -= Utils::calcDist(res.output.route[a].pos, res.output.route[a-1].pos);
        res.length -= Utils::calcDist(res.output.route[a].pos, res.output.route[a+1].pos);
        res.length -= Utils::calcDist(res.output.route[b].pos, res.output.route[b-1].pos);
        res.length -= Utils::calcDist(res.output.route[b].pos, res.output.route[b+1].pos);
        swap(res.output.route[a], res.output.route[b]);
        assert(res.output.route_id[0][va.id] < res.output.route_id[1][va.id]);
        assert(res.output.route_id[0][vb.id] < res.output.route_id[1][vb.id]);
        res.length += Utils::calcDist(res.output.route[a].pos, res.output.route[a-1].pos);
        res.length += Utils::calcDist(res.output.route[a].pos, res.output.route[a+1].pos);
        res.length += Utils::calcDist(res.output.route[b].pos, res.output.route[b-1].pos);
        res.length += Utils::calcDist(res.output.route[b].pos, res.output.route[b+1].pos);
        res.score = Utils::calcScore(res.length);
    }
    return res;
}

State State::generateState(const State& input_state) {
    State res = input_state;
    int x = ryuka.rand(res.output.nodes.size());
    int y = ryuka.rand(input.n);
    int count = 0;
    for(int i : res.output.nodes) count += (i == y); 
    if(count == 0) {
        res.output.nodes[x] = y;
        return initState(res.output.nodes);
    }
    return res;
}

State2 State2::convert(const State& input_state) {
    State2 res;
    res.output = input_state.output;
    res.score = input_state.score;
    res.length = input_state.length;
    return res;
}

int Utils::calcDist(const pair<int,int>& a, const pair<int,int>& b) {
    int res = abs(a.first - b.first) + abs(a.second - b.second);
    return res;
}

/*TODO: ここでスコアを計算する*/
long long Utils::calcLength(const Output& output) {
    long long res = 0;
    for(int i = 0; i < output.route.size() - 1; i++) {
        res += Utils::calcDist(output.route[i].pos, output.route[i+1].pos);
    }
    return res;
}

long long Utils::calcScore(long long length) {
    long long res = (100000000 / (1000 + length));
    return res;
}

long long Utils::calcScore(const Output& output) {
    long long res = (100000000 / (1000 + Utils::calcLength(output)));
    return res;
}

Output Utils::runInsertTSP(const vector<int>& nodes) {
    vector<Node> route;
    vector<vector<int>> route_id(2, vector<int>(input.n));
    map<pair<int,int>,int> node_counter;
    route.push_back(Node(make_pair(400, 400), 0, -1));
    route.push_back(Node(make_pair(400, 400), 1, -1));
    for(int id : nodes) {
        vector<pair<int,int>> positions = {input.src(id), input.dst(id)};
        int src_pos = 0;
        for(int p = 0; p < 2; p++) {
            pair<int,int> pos = positions[p];
            int min_dist = 1<<30, min_id = -1;
            for(int j = 0; j < route.size()-1; j++) {
                int dist = Utils::calcDist(route[j].pos, pos) + Utils::calcDist(route[j+1].pos, pos);
                if(dist < min_dist && j >= src_pos) {
                    min_dist = dist;
                    min_id = j;
                }
            }
            assert(min_id != -1);
            route.insert(route.begin() + min_id + 1, Node(pos, p, id));
            src_pos = min_id + 1;
        }
    }
    for(int i = 1; i < route.size()-1; i++) {
        Node v = route[i];
        if(v.id >= 0) route_id[v.p][v.id] = i;
    }
    Output res;
    res.nodes = nodes;
    res.route = route;
    res.route_id = route_id;
    return res;
}

int main(int argc, char* argv[]) {
    toki.init();
    input.read();   
    IterationControl<State> sera;
    State pre_ans = sera.climb(1.4, State::initState());
    IterationControl<State2> funaq;
    State2 ans = funaq.climb(1.8 - toki.elapsed(), State2::convert(pre_ans));
    //State ans = State::initState();
    ans.output.print();
    cerr << "[INFO] - main - MyScore = " << ans.score << "\n";
    cerr << "[INFO] - main - MyLength = " << ans.length << "\n";
}
