﻿//------------------------------------------------------------------------------
/// @file
/// @author   ハル研究所プログラミングコンテスト実行委員会
///
/// @copyright  (C)HAL Laboratory, Inc.
/// @attention  このファイルの利用は、同梱のREADMEにある
///             利用条件に従ってください。
//------------------------------------------------------------------------------

#include "Answer.hpp"
#include <map>
#include <unordered_map>
#include <queue>
#include <cmath>
#include <algorithm>
#include <tuple>

template<typename T>
void hash_combine(size_t & seed, T const& v) {
    //基本型に関するハッシュ生成は標準ライブラリが提供している
    std::hash<T> primitive_type_hash;

    //生成したハッシュを合成する。このコードはboostものを使用する
    seed ^= primitive_type_hash(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace hpc {
bool operator==(const Vector2& lhs, const Vector2& rhs){
    return lhs.x == rhs.x && lhs.y == rhs.y;
}
}

namespace std {
template<>
class hash<hpc::Vector2>{
    public:
        size_t operator()(const hpc::Vector2& data) const {
            std::size_t seed = 0;
            hash_combine(seed, data.x);
            hash_combine(seed, data.y);

            return seed;
        }
};
}

//------------------------------------------------------------------------------
namespace hpc {

typedef std::vector<Vector2> Ver;
typedef std::vector<Vector2> Path;
// {対象までの距離, 座標}
typedef std::pair<double, Vector2> P_dist;
typedef std::pair<int, int> P;

const float EPS = 1e-5;
//const double EPS = 1e-3;
const double INF = 1e50;
/// rabbit : 0, others : 1 - 20
const int MAX_V = 1 + Parameter::MaxScrollCount;

/// ウサギ&巻物 間の速さをジャンプ力を1.0とした場合の移動時間
double distance[25][25];

/// dijkstra方で用いる距離
int dist_scroll[Parameter::MaxScrollCount][Parameter::StageHeight][Parameter::StageWidth];

/// 巻物の取る順番
std::vector<int> targets;

struct Edge{int y, x, weight;};

/// 各タイルを結んだグラフ {nx, ny, weight}
std::vector<Edge> Graph[Parameter::StageHeight][Parameter::StageWidth];

double dist(Vector2 v1, Vector2 v2);
Vector2 calc_point_rot(Vector2 v1, Vector2 v2, double theta);
double dist_opt(const Stage& aStage, Vector2 v1, Vector2 v2, int top_k, const std::vector<double>& Theta);
void distance_init(const Stage& aStage);
void warshall_floyd(Ver vertices);
void calc_distance_greedy(Ver vertices);
void calc_distance_dijkstra(const Stage& aStage, Ver vertices);
void calc_distance_search_topk(const Stage& aStage, Ver vertices, int topk, const std::vector<double>& Theta);
void sales_man_init(int size);
double sales_man(int S, int v, int size);
void vertices_init(Ver& vertices, const Stage& aStage);
void build_target_sequence();
Path get_path_opt(const Stage& aStage, Vector2 v1, Vector2 v2, int top_k, const std::vector<double>& Theta);
void make_graph(const Stage& aStage);
bool same_float(Vector2 v1, Vector2 v2);
bool same(Vector2 v1, Vector2 v2);


////////////////////////////////////////////////////////////////////////

/// v1, v2間の直線距離
double dist(Vector2 v1, Vector2 v2){
    double diff_y = v1.y - v2.y;
    double diff_x = v1.x - v2.x;
    return std::sqrt(diff_y*diff_y + diff_x*diff_x);
}

/// v1からv2へrot回転された点を算出
Vector2 calc_point_rot(Vector2 v1, Vector2 v2, double theta){
    Vector2 target = {0, 0};
    float vec_x = v2.x - v1.x;
    float vec_y = v2.y - v1.y;
    target.x = v1.x + vec_x * cos(theta) - vec_y * sin(theta);
    target.y = v1.y + vec_x * sin(theta) + vec_y * cos(theta);
    if(target.y < 0){
        target.x = v1.x + (target.x - v1.x) * (v1.y / (v1.y - target.y));
        target.y = 0;
    }
    if(target.x < 0){
        target.y = v1.y + (target.y - v1.y) * (v1.x / (v1.x - target.x));
        target.x = 0;
    }
    if(target.y >= Parameter::StageWidth){
        const float Y = Parameter::StageWidth - EPS;
        target.x = v1.x + (target.x - v1.x) * ((v1.y - Y) / (v1.y - target.y));
        target.y = Y;
    }
    if(target.x >= Parameter::StageHeight){
        const float X = Parameter::StageHeight - EPS;
        target.y = v1.y + (target.y - v1.y) * ((v1.x - X) / (v1.x - target.x));
        target.x = X;
    }
    return target;
}

/// v1とv2が同じ位置かどうか判定
bool same_float(Vector2 v1, Vector2 v2){
    return (fabs(v1.x - v2.x) < EPS && fabs(v1.y - v2.y) < EPS);
}

/// v1とv2が同じマスかどうか判定
bool same(Vector2 v1, Vector2 v2){
    return ((int)v1.x == (int)v2.x && (int)v1.y == (int)v2.y);
}

// TODO : check efficiency
template<typename T>
void queue_clear(T &que)
{
   T empty;
   std::swap(que, empty);
}

class Comp_P_dist{
    public:
    bool operator() (const P_dist& a, const P_dist& b){
        return a.first > b.first;
    }
};

/// v1, v2間の移動にかかるターン数(ジャンプ力 = 1)
double dist_opt(const Stage& aStage, Vector2 v1, Vector2 v2, int top_k, const std::vector<double>& Theta){
    std::priority_queue<P_dist, std::vector<P_dist>, Comp_P_dist> que;
    que.push(P_dist(dist(v1, v2), v1));
    //const std::vector<double> Theta = {M_PI/2, M_PI/3, M_PI/6, 0, -M_PI/6, -M_PI/3, -M_PI/2};
    int turn = 0;
    while(1){
        Ver pos_array;
        bool finish = false;
        while(!que.empty() && (int)pos_array.size() < top_k){
            auto x = que.top(); que.pop();
            pos_array.push_back(x.second);
            //if(x.first < EPS){
            auto pos = x.second;
            if((int)pos.x == (int)v2.x && (int)pos.y == (int)v2.y){
                //printf("Found : turn = %d\n", turn);
                finish = true;
                break;
            }
        }

        if(finish)
            break;

        turn++;
        queue_clear<decltype(que)>(que);
        for(auto v : pos_array){
            for(auto theta : Theta){
                // +theta, -thetaの両方を行う
                int nloop = (theta == 0) ? 1 : 2;
                for(int i = 0; i < nloop; i++){
                    theta = (i == 0) ? theta : -theta;
                    auto target = calc_point_rot(v, v2, theta);
                    auto nv = aStage.getNextPos(v, 1.0, target);
                    que.push(P_dist(dist(nv, v2), nv));
                }
            }
        }
    }
    return turn;
}

class Vector2Compare{
    public:
        bool operator()(const Vector2& a, const Vector2& b){return (a.x == b.x) ? a.y < b.y : a.x < b.x;};
};

/// v1, v2間の移動にかかるターン数が最小となるパス
Path get_path_opt(const Stage& aStage, Vector2 v1, Vector2 v2, int top_k, const std::vector<double>& Theta){
    std::priority_queue<P_dist, std::vector<P_dist>, Comp_P_dist> que;
    que.push(P_dist(dist(v1, v2), v1));
    // 前回の座標値を記録
    //std::map<Vector2, Vector2, Vector2Compare> pre_pos;
    //std::unordered_map<MyVector2, Vector2, MyVector2::Hash> pre_pos;
    std::unordered_map<Vector2, Vector2> pre_pos;
    Vector2 last;
    while(1){
        Ver pos_array;
        bool finish = false;
        while(!que.empty() && (int)pos_array.size() < top_k){
            auto x = que.top(); que.pop();
            pos_array.push_back(x.second);
            //TODO : check validity
            //if(x.first < EPS){
            auto pos = x.second;
            if((int)pos.x == (int)v2.x && (int)pos.y == (int)v2.y){
                //printf("Found : turn = %d\n", turn);
                last = x.second;
                finish = true;
                break;
            }
        }

        if(finish)
            break;

        //TODO : need this??
        queue_clear<decltype(que)>(que);
        for(auto v : pos_array){
            for(auto theta : Theta){
                // +theta, -thetaの両方を行う
                int nloop = (theta == 0) ? 1 : 2;
                for(int i = 0; i < nloop; i++){
                    theta = (i == 0) ? theta : -theta;
                    auto target = calc_point_rot(v, v2, theta);
                    auto nv = aStage.getNextPos(v, aStage.rabbit().power(), target);
                    if(pre_pos.find(nv) != pre_pos.end()){
                        //printf("Found!!\n");
                        continue;
                    }
                    que.push(P_dist(dist(nv, v2), nv));
                    pre_pos[nv] = v;
                }
            }
        }
    }
    Path path;
    while(!(last.x == v1.x && last.y == v1.y)){
        path.push_back(last);
        last = pre_pos[last];
    }
    reverse(path.begin(), path.end());
    return path;
}

/// 各点間をINFに初期化
void distance_init(const Stage& aStage){
    int size = aStage.scrolls().count();
    for(int i = 0; i < size; i++)
        for(int j = 0; j < size; j++)
            if(i == j) distance[i][j] = 0;
            else distance[i][j] = INF;
}

/// 全点間距離を求める。()
void warshall_floyd(Ver vertices){
    int size = vertices.size();
    for(int i = 0; i < size; i++)
        for(int j = 0; j < size; j++)
            if(i == j) distance[i][j] = 0;
            else distance[i][j] = INF;

    for(int k = 0; k < size; k++)
        for(int i = 0; i < size; i++)
            for(int j = 0; j < size; j++)
                distance[i][j] = std::min(distance[i][j], distance[i][k] + distance[k][j]);
}

/// 直線距離で全点間のdistanceを更新
void calc_distance_greedy(Ver vertices){
    int size = vertices.size();
    for(int i = 0; i < size; i++){
        for(int j = i+1; j < size; j++){
            auto v1 = vertices[i];
            auto v2 = vertices[j];
            distance[i][j] = distance[j][i] = dist(v1, v2);
        }
    }
}

/// ビームサーチによる経路探索
void calc_distance_search_topk(const Stage& aStage, Ver vertices, int topk, const std::vector<double>& Theta){
    int size = vertices.size();
    for(int i = 0; i < size; i++){
        for(int j = i+1; j < size; j++){
            auto v1 = vertices[i];
            auto v2 = vertices[j];
            distance[i][j] = distance[j][i] = dist_opt(aStage, v1, v2, topk, Theta);
        }
    }
}

/// dijkstra法を基にdistanceを更新
void calc_distance_dijkstra(const Stage& aStage, Ver vertices){
    int size = vertices.size();
    // i = 0 : rabbit
    for(int i = 1; i < size; i++){
        for(int j = 0; j < size; j++){
            //auto v1 = vertices[i];
            int idx = i-1;
            auto v2 = vertices[j];
            distance[i][j] = distance[j][i] = dist_scroll[idx][(int)v2.y][(int)v2.x];
        }
    }
}

/// 最短経路での移動先の点
int next_vertex[1 << MAX_V][MAX_V];
/// 計算に使用 TODO:change
double DP[1 << MAX_V][MAX_V];
/// DP配列の初期化
void sales_man_init(int size){
    for(int i = 0; i < (1 << size); i++){
        for(int j = 0; j < size; j++){
            DP[i][j] = -1.0;
            next_vertex[i][j] = -1;
        }
    }
}

/// すべての頂点を一度ずつめぐり帰ってくる経路のうち、重みの総和の最小値を求める。
/// O(2**n n**2)
double sales_man(int S, int v, int size){
    if(DP[S][v] >= 0){
        return DP[S][v];
    }

    //if(S == (1 << size)-1 && v == 0)
    if(S == (1 << size)-1)
        return DP[S][v] = 0;

    double res = INF;
    for(int u = 0; u < size; u++){
        if(!(S >> u & 1)){
            double tmp = sales_man(S | 1 << u, u, size) + distance[v][u];
            if(res > tmp){
                res = tmp;
                next_vertex[S][v] = u;
            }
        }
    }
    return DP[S][v] = res;
}

/// すべての頂点を一度ずつめぐり帰ってくる経路のうち、重みの総和の最小値を求める。
/// TODO : 巻物の取得数も考慮する.  ほぼ効果なし。
/// O(2**n n**2)
double sales_man_scroll(int S, int v, int size){
    if(DP[S][v] >= 0){
        return DP[S][v];
    }

    //if(S == (1 << size)-1 && v == 0)
    if(S == (1 << size)-1)
        return DP[S][v] = 0;

    int pop_count = 0;
    for(int i = 0; i < size; i++)
        if((S >> i) & 1)
            pop_count++;
    float decay = 1.0;
    for(int p = 0; p < pop_count-1; p++)
        decay *= 1.1;

    double res = INF;
    for(int u = 0; u < size; u++){
        if(!(S >> u & 1)){
            double tmp = sales_man(S | 1 << u, u, size) + (float)distance[v][u] / decay;
            if(res > tmp){
                res = tmp;
                next_vertex[S][v] = u;
            }
        }
    }
    return DP[S][v] = res;
}

/// vertice配列の初期化
void vertices_init(Ver& vertices, const Stage& aStage){
    vertices.push_back(aStage.rabbit().pos());
    for(auto scroll : aStage.scrolls())
        vertices.push_back(scroll.pos());
}

/// next_vertex配列からtargets配列を生成
void build_target_sequence(){
    int v = 0;
    int S = 0;
    while(1){
        S |= 1 << v;
        v = next_vertex[S][v];
        if(v == -1) break;
        /// v=0はRabbitであるため。
        targets.push_back(v-1);
    }

    //printf("%d : ", (int)targets.size());
    //for(auto idx : targets){
    //    printf("%d ", idx);
    //}
    //printf("\n");
}

/// 地形情報からGraphを初期化.(いらないかも)
void make_graph(const Stage& aStage){
    for(int y = 0; y < Parameter::StageHeight; y++){
        for(int x = 0; x < Parameter::StageWidth; x++){
            Graph[y][x].clear();
            auto terrain = aStage.terrain(Vector2{(float)x, (float)y});
            int weight = 100000;
            if(terrain == Terrain::Plain)
                weight = 3;
            else if(terrain == Terrain::Bush)
                weight = 5;
            else if(terrain == Terrain::Sand)
                weight = 10;
            else if(terrain == Terrain::Pond)
                weight = 30;

            const int dx[] = {-1, 1, 0, 0};
            const int dy[] = {0, 0, 1, -1};
            for(int i = 0; i < 4; i++){
                int ny = y + dy[i];
                int nx = x + dx[i];
                if(ny < 0 || nx < 0) continue;
                if(ny >= Parameter::StageHeight) continue;
                if(nx >= Parameter::StageWidth) continue;
                Graph[y][x].push_back(Edge{ny, nx, weight});
            }
        }
    }
}

///ダイクストラ法
void dijkstra(const Stage& aStage, int s){
    const int Y = Parameter::StageHeight;
    const int X = Parameter::StageWidth;
    typedef std::tuple<int, int, int> T;
    for(int y = 0; y < Y; y++)
        for(int x = 0; x < X; x++)
            dist_scroll[s][y][x] = 1e6;
    auto scrolls = aStage.scrolls(); 
    int y1 = (int)scrolls[s].pos().y, x1 = (int)scrolls[s].pos().x;
    dist_scroll[s][y1][x1] = 0;
    // {dist, y, x}
    std::priority_queue<T, std::vector<T>, std::greater<T>> que;
    que.push(T(0, y1, x1));

    while(!que.empty()){
        auto t = que.top(); que.pop();
        int d = std::get<0>(t);
        int y = std::get<1>(t);
        int x = std::get<2>(t);
        if(dist_scroll[s][y][x] < d) continue;
        for(auto&& edge : Graph[y][x]){
            int ny = edge.y;
            int nx = edge.x;
            int nd = d + edge.weight;
            if(dist_scroll[s][ny][nx] > nd){
                dist_scroll[s][ny][nx] = nd;
                que.push(T(nd, ny, nx));
            }
        }
    }
}

/// srcから見たdestが第何象限にあるか
int get_direction(Vector2 src, Vector2 dest){
    float x1 = src.x;
    float y1 = src.y;
    float x2 = dest.x;
    float y2 = dest.y;

    float vec_x = x2 - x1;
    float vec_y = y2 - y1;

    if(vec_x >= 0 && vec_y >= 0)
        return 0;
    else if(vec_x < 0 && vec_y >= 0)
        return 1;
    else if(vec_x < 0 && vec_y < 0)
        return 2;
    else if(vec_x >= 0 && vec_y < 0)
        return 3;
    else
        return -1;
}

/// {y1, x1} から scrolls[dest] までの経路を復元.
/// TODO:最後はdest2に近い方角で終わる。 ラストから決めていかないとまずい？？
void get_path_from_dijkstra(const Stage& aStage, Path& path, int y1, int x1, int dest, int dest2){
    //TODO : 最初の地点を保存？
    int d = dist_scroll[dest][y1][x1];
    auto target = aStage.scrolls()[dest].pos();
    int y = y1, x = x1;
    while(d != 0){
        //TODO : to center? 
        // 特にtargetが池の場合は少し入ってすぐ出るべきだから中心はまずい
        //path.push_back(Vector2{(float)x+(float)0.5, (float)y+(float)0.5});
        int dx[2][4][4] = {{{1, 0, -1, 0}, {0, -1, 0, 1}, {-1, 0, 1, 0}, {0, 1, 0, -1}},
                           {{0, -1, 0, 1}, {-1, 0, 1, 0}, {0, 1, 0, -1}, {1, 0, -1, 0}}};
        int dy[2][4][4] = {{{0, 1, 0, -1}, {1, 0, -1, 0}, {0, -1, 0, 1}, {-1, 0, 1, 0}},
                           {{1, 0, -1, 0}, {0, -1, 0, 1}, {-1, 0, 1, 0}, {0, 1, 0, -1}}};
        int direction = get_direction(Vector2{(float)x+(float)0.5, (float)y+(float)0.5}, target);
        const int turn = path.size();

        for(int i = 0; i < 4; i++){
            //int ny = y + dy[i];
            //int nx = x + dx[i];
            int ny = y + dy[turn%2][direction][i];
            int nx = x + dx[turn%2][direction][i];
            if(ny < 0 || nx < 0) continue;
            if(ny >= Parameter::StageHeight) continue;
            if(nx >= Parameter::StageWidth) continue;
            //Graph[y][x].push_back(Edge{ny, nx, weight});
            int nd = dist_scroll[dest][ny][nx];

            auto terrain = aStage.terrain(Vector2{(float)nx, (float)ny});
            int weight = 100000;
            if(terrain == Terrain::Plain)
                weight = 3;
            else if(terrain == Terrain::Bush)
                weight = 5;
            else if(terrain == Terrain::Sand)
                weight = 10;
            else if(terrain == Terrain::Pond)
                weight = 30;
            if(nd + weight == d){
                float tx = x + 0.5, ty = y + 0.5;
                if(nd == 0){
                    tx = ((1+EPS)*(nx+0.5) + (x+0.5)) / (2. + EPS);
                    ty = ((1+EPS)*(ny+0.5) + (y+0.5)) / (2. + EPS);
                }else{
                    auto ter_org = aStage.terrain(Vector2{x+EPS, y+EPS});
                    auto ter_nxt = aStage.terrain(Vector2{nx+EPS, ny+EPS});
                    if(ter_org == ter_nxt){
                        tx = (x + nx + 1) / 2.0;
                        ty = (y + ny + 1) / 2.0;
                    }
                    else if(ter_org < ter_nxt){
                        tx = ((nx+0.5) + (1+EPS)*(x+0.5)) / (2. + EPS);
                        ty = ((ny+0.5) + (1+EPS)*(y+0.5)) / (2. + EPS);
                    }
                    else{
                        tx = ((1+EPS)*(nx+0.5) + (x+0.5)) / (2. + EPS);
                        ty = ((1+EPS)*(ny+0.5) + (y+0.5)) / (2. + EPS);
                    }
                }
                path.push_back(Vector2{tx, ty});
                d = nd, y = ny, x = nx;
                break;
            }
        }
    }
}


//------------------------------------------------------------------------------
/// コンストラクタ
/// @detail 最初のステージ開始前に実行したい処理があればここに書きます
Answer::Answer()
{
}

//------------------------------------------------------------------------------
/// デストラクタ
/// @detail 最後のステージ終了後に実行したい処理があればここに書きます
Answer::~Answer()
{
}

// for getTragetPos
Path path;
int path_idx = 0;
bool big_jump;

/// TODO : Game Emulation
/// targets配列の順序をランダムに変えてみる？

//------------------------------------------------------------------------------
/// 各ステージ開始時に呼び出される処理
/// @detail 各ステージに対する初期化処理が必要ならここに書きます
/// @param aStage 現在のステージ
//int test;
void Answer::initialize(const Stage& aStage)
{
    //printf("test : %d\n", test++);
    //if(aStage.turn() < 1)
    //    return;

    // {rabbit, scroll0, scroll1, ...}
    Ver vertices;
    vertices_init(vertices, aStage);

    make_graph(aStage);
    for(int v = 0; v < (int)aStage.scrolls().count(); v++)
        dijkstra(aStage, v);
    calc_distance_dijkstra(aStage, vertices);

    //calc_distance_greedy(vertices);
    //const std::vector<double> Theta = {M_PI*sqrt(3.0/4), M_PI/6, 0};
    //const int topk = 10;
    //calc_distance_search_topk(aStage, vertices, topk, Theta);

    sales_man_init(vertices.size());
    //sales_man(1 << 0, 0, vertices.size());
    sales_man_scroll(1 << 0, 0, vertices.size());

    build_target_sequence();

    path.clear();
    big_jump = false;
}

//------------------------------------------------------------------------------
/// 毎フレーム呼び出される処理
/// @detail 移動先を決定して返します
/// @param aStage 現在のステージ
/// @return 移動の目標座標
Vector2 Answer::getTargetPos(const Stage& aStage)
{
    auto pos = aStage.rabbit().pos();
/*
    for(auto scroll : aStage.scrolls()) {
        // まだ手に入れていない巻物を探して、そこに向かって飛ぶ
        if (!scroll.isGotten()) {
            return scroll.pos();
        }
    }
    return pos;
*/
    if(aStage.turn() == 0 || same(pos, path[path.size()-1])){
        auto scrolls = aStage.scrolls();
        for(int i = 0; i < (int)targets.size(); i++){
            int idx = targets[i];
            auto scroll = scrolls[idx];
            if (!scroll.isGotten()) {
                path.clear();
                if(i == (int)targets.size()-1)
                    get_path_from_dijkstra(aStage, path, pos.y, pos.x, idx, -1);
                else
                    get_path_from_dijkstra(aStage, path, pos.y, pos.x, idx, targets[i+1]);
                path_idx = 0;
                break;
            }
        }
    }

    Vector2 target;
    Terrain tar_terrain;
    Vector2 pre_point = pos;
    Terrain pre_terrain = aStage.terrain(pre_point);
    bool first = true;
    int tmp_idx = 0;

    for(int k = std::min(path_idx + 5, (int)path.size()-1); k > path_idx; k--){
        auto tmp_target = path[k];
        if(same(tmp_target, pos)){
            path_idx = k;
            break;
        }
    }

    while(1){
        // TODO : 逆戻り解消
        target = path[path_idx + tmp_idx];
        tar_terrain = aStage.terrain(target);

        auto next = aStage.getNextPos(pos, aStage.rabbit().power(), target);
        auto nex_terrain = aStage.terrain(next);

        if(!first && tar_terrain > pre_terrain){
            target = pre_point;
            break;
        }
        else if(path_idx+1 < (int)path.size() && same_float(target, next)){
            path_idx++;
        }
        else if(path_idx+1 < (int)path.size() && same(target, next)){
            tmp_idx++;
        }
        else{
            if(nex_terrain == tar_terrain)
                target = next;
            else if(!first)
                target = pre_point;
            break;
        }

        pre_point = target;
        pre_terrain = tar_terrain;
        first = false;
    }
    return target;

/* //池ポチャ
    const Vector2 fail = {-1, -1};
    Vector2 target = fail;
    Terrain tar_terrain;
    Vector2 pre_point = pos;
    Terrain pre_terrain = aStage.terrain(pre_point);
    //bool selected = false;

    for(int i = (int)path.size()-1; i >= path_idx; i--){
        auto tmp = path[i];
        if(same(pos, tmp)){
            path_idx = i;
            break;
        }
    }

    if(same_float(pos, path[path_idx]))
        path_idx++;

    for(int k = 0; k < 5; k++){
        if(path_idx + k >= (int)path.size()){
            //selected = true;
            break;
        }
        target = path[path_idx + k];
        tar_terrain = aStage.terrain(target);
        if((k > 0 && tar_terrain > pre_terrain){
            target = pre_point;
            //selected = true;
            break;
        }
        pre_point = target;
        pre_terrain = tar_terrain;
    }

    if(same_float(fail, target))
        target = path[path.size()-1];
    return target;
*/
}

//------------------------------------------------------------------------------
/// 各ステージ終了時に呼び出される処理
/// @detail 各ステージに対する終了処理が必要ならここに書きます
/// @param aStage 現在のステージ
void Answer::finalize(const Stage& aStage)
{
    targets.clear();
}

} // namespace
// EOF
