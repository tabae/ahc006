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

struct Output {
    /*TODO: ここに出力変数を定義する*/
    vector<int> nodes;
    vector<pair<int,int>> root;
    vector<vector<int>> root_id;
    Output();
    void print();
};

/*解を管理するクラス*/
struct State {
    Output output;
    long long score;
    State() : score(0) {}
    static State initState();
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
    long long calcScore(const Output& output);
    int calcDist(const pair<int,int>& a, const pair<int,int>& b);
    pair<vector<pair<int,int>>,vector<vector<int>>> runInsertTSP(const vector<int>& nodes);
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
    cout << root.size();
    for(pair<int,int> node: root) cout << " " << node.first << " " << node.second;
    cout << endl;
}

/*TODO: ここで初期解を作成する*/
State State::initState() {
    State res;
    // ランダムに50個選ぶ
    vector<bool> seen(input.n+1, false);
    while(res.output.nodes.size() < 50) {
        int id = ryuka.rand(input.n);
        if(seen[id]) continue;
        res.output.nodes.push_back(id);
        seen[id] = true;
    }
    // とりあえず、挿入法で適当な巡回経路を作る。
    auto [root, root_id] = Utils::runInsertTSP(res.output.nodes);
    res.output.root = root;
    res.output.root_id = root_id;
    res.score = Utils::calcScore(res.output);
    return res;
}

/*TODO: ここでinput_stateを変化させた解を作る（局所探索）*/
State State::generateState(const State& input_state) {
    State res = input_state;
    res.score = Utils::calcScore(res.output);
    return res;
}

State2 State2::generateState(const State2& input_state) {
    State2 res = input_state;
    res.score = Utils::calcScore(res.output);
    return res;
}

State2 State2::convert(const State& input_state) {
    State2 res;
    res.output = input_state.output;
    res.score = input_state.score;
    return res;
}

int Utils::calcDist(const pair<int,int>& a, const pair<int,int>& b) {
    int res = abs(a.first - b.first) + abs(a.second - b.second);
    return res;
}


/*TODO: ここでスコアを計算する*/
long long Utils::calcScore(const Output& output) {
    long long res = 0;
    for(int i = 0; i < output.root.size() - 1; i++) {
        res += Utils::calcDist(output.root[i], output.root[i+1]);
    }
    return res;
}

pair<vector<pair<int,int>>, vector<vector<int>>> Utils::runInsertTSP(const vector<int>& nodes) {
    vector<pair<int,int>> root;
    vector<vector<int>> root_id(2, vector<int>(nodes.size()));
    root.push_back({400, 400});
    root.push_back({400, 400});
    for(int i = 0; i < nodes.size(); i++) {
        vector<pair<int,int>> positions = {input.src(nodes[i]), input.dst(nodes[i])};
        for(int p = 0; p < 2; p++) {
            pair<int,int> pos = positions[p];
            int min_dist = 1<<30, min_id = -1;
            for(int j = 0; j < root.size()-1; j++) {
                int dist = Utils::calcDist(root[j], pos) + Utils::calcDist(root[j+1], pos);
                if(dist < min_dist) {
                    if(p == 1 && j < root_id[0][i]) ;
                    else {
                        min_dist = dist;
                        min_id = j;
                    }
                }
            }
            assert(min_id != -1);
            root.insert(root.begin() + min_id + 1, pos);
            root_id[p][i] = min_id + 1;
        }
    }
    return {root, root_id};
 }

int main(int argc, char* argv[]) {
    toki.init();
    input.read();   
    IterationControl<State> sera;
    //State ans = sera.climb(1.8, State::initState());
    State ans = State::initState();
    ans.output.print();
    cerr << "[INFO] - main - MyScore = " << ans.score << "\n";
}
