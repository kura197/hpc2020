//------------------------------------------------------------------------------
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
#include <cassert>
#include <random>

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

float POWER[30];

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
void build_target_sequence_from_tsp();
void make_graph(const Stage& aStage);
bool same_float(Vector2 v1, Vector2 v2);
bool same(Vector2 v1, Vector2 v2);
int get_terrain_weight(const Stage& aStage, Vector2 pos);

unsigned int randxor();
void get_path_from_dijkstra(const Stage& aStage, Path& path, int y1, int x1, int src, int dest, int dest2);
/// ;y1, x1} から scrolls[dest] までの経路を復元. 実装を簡単にするため、目標値はすべてマスの中心にする
void get_path_from_dijkstra_simple(const Stage& aStage, Path& path, int y1, int x1, int src, int dest, int dest2);
/// ウサギが次にジャンプすべき座標を計算
Vector2 next_path(const Stage& aStage, Path path, int& path_idx, Vector2 cur_pos, float power);
/// dijkstra_pathを全点間で作成
void build_all_dijkstra_path(const Stage& aStage);
// v1 -> v2 へのウサギの通る座標(sequence_simple)を更新
void update_path_simple(const Stage& aStage, int v1, int v2, float power);
/// targets配列の初期化 & 
void targets_init(const Stage& aStage);
/// targetsからpathをupdate
void update_path_from_targets(const Stage& aStage);
/// 終了までに要するターン数を計算
int path_length(const Stage& aStage);
/// 終了までに要するターン数を計算
float path_length_from_dijkstra(const Stage& aStage);
/// targetsのv1番目とv2番目を交換し、再計算
void change_vertex(const Stage& aStage, int v1, int v2);
/// targetsのv1番目とv2番目を交換し、再計算
void change_vertex_simple(const Stage& aStage, int v1, int v2);
/// 2-opt法でTSPを解く (1.404sec)
void tsp_2opt(const Stage& aStage);
/// 2-opt法でTSPを解く (1.404sec)
void tsp_2opt_simple(const Stage& aStage);
Vector2 MygetTargetPos(const Stage& aStage);


////////////////////////////////////////////////////////////////////////

class EmulateGame{
    private:
        Stage stage;

    public:
        EmulateGame(const Stage& aStage){
            stage = aStage;
        }

        int run(){
            int current_turn = 0;
            while (!stage.isEnd() && stage.turn() < Parameter::GameTurnLimit && !stage.isOutOfBounds(stage.rabbit().pos())) {
                // ターン開始
                auto targetPos = MygetTargetPos(stage);
                stage.update(targetPos);
                stage.advanceTurn();
                current_turn++;
            }
            return current_turn;
        }
};


int get_terrain_weight(const Stage& aStage, Vector2 pos){
    int weight = 0;
    auto terrain = aStage.terrain(pos);
    if(terrain == Terrain::Plain)
        //weight = 3;
        weight = 6;
    else if(terrain == Terrain::Bush)
        //weight = 5;
        weight = 10;
    else if(terrain == Terrain::Sand)
        //weight = 10;
        weight = 20;
    else if(terrain == Terrain::Pond)
        //weight = 30;
        //weight = 60;
        weight = 60;
    return weight;
}

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
void build_target_sequence_from_tsp(){
    int v = 0;
    int S = 0;
    targets.clear();
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
            int weight = get_terrain_weight(aStage, Vector2{(float)x, (float)y}); 

            const int dx[] = {-1, 1, 0, 0};
            const int dy[] = {0, 0, 1, -1};
            for(int i = 0; i < 4; i++){
                int ny = y + dy[i];
                int nx = x + dx[i];
                if(ny < 0 || nx < 0) continue;
                if(ny >= Parameter::StageHeight) continue;
                if(nx >= Parameter::StageWidth) continue;

                // TODO : expriment
                int new_weight = weight + get_terrain_weight(aStage, Vector2{(float)nx, (float)ny});
                
                //Graph[y][x].push_back(Edge{ny, nx, weight});
                Graph[y][x].push_back(Edge{ny, nx, new_weight});
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
Vector2 first_step[MAX_V];
void get_path_from_dijkstra(const Stage& aStage, Path& path, int y1, int x1, int src, int dest, int dest2){
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

            int weight = get_terrain_weight(aStage, Vector2{(float)nx, (float)ny});

            //TODO : experiment
            weight += get_terrain_weight(aStage, Vector2{(float)x, (float)y});

            if(nd + weight == d){
                float tx = x + 0.5, ty = y + 0.5;
                if(nd == 0){
                    if(dest2 < 0){
                        tx = ((1+EPS)*(nx+0.5) + (x+0.5)) / (2. + EPS);
                        ty = ((1+EPS)*(ny+0.5) + (y+0.5)) / (2. + EPS);
                    }
                    else{
                        auto next_step = first_step[dest2];
                        tx = ((1+EPS)*(nx+0.5) + (next_step.x+0.5)) / (2. + EPS);
                        ty = ((1+EPS)*(ny+0.5) + (next_step.y+0.5)) / (2. + EPS);
                    }
                }else{
                    auto ter_org = aStage.terrain(Vector2{x+EPS, y+EPS});
                    auto ter_nxt = aStage.terrain(Vector2{nx+EPS, ny+EPS});
                    if(ter_org == ter_nxt){
                        tx = (x + nx + 1) / 2.0;
                        ty = (y + ny + 1) / 2.0;
                    }
                    else if(ter_org < ter_nxt){
                        ////TODO : もう少し簡単に
                        tx = ((nx+0.5) + (1+EPS)*(x+0.5)) / (2. + EPS);
                        ty = ((ny+0.5) + (1+EPS)*(y+0.5)) / (2. + EPS);
                    }
                    else{
                        tx = ((1+EPS)*(nx+0.5) + (x+0.5)) / (2. + EPS);
                        ty = ((1+EPS)*(ny+0.5) + (y+0.5)) / (2. + EPS);
                    }

                    if(path.size() == 0 && src >= 0)
                        first_step[src] = Vector2{tx, ty};
                }
                path.push_back(Vector2{tx, ty});
                d = nd, y = ny, x = nx;
                break;
            }
        }
    }
}

/// {y1, x1} から scrolls[dest] までの経路を復元. 実装を簡単にするため、目標値はすべてマスの中心にする
void get_path_from_dijkstra_simple(const Stage& aStage, Path& path, int y1, int x1, int src, int dest, int dest2){
    int d = dist_scroll[dest][y1][x1];
    auto target = aStage.scrolls()[dest].pos();
    int y = y1, x = x1;
    while(d != 0){
        //TODO : to center? 
        // 特にtargetが池の場合は少し入ってすぐ出るべきだから中心はまずい
        path.push_back(Vector2{(float)x+(float)0.5, (float)y+(float)0.5});
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

            int weight = get_terrain_weight(aStage, Vector2{(float)nx, (float)ny});

            //TODO : experiment
            weight += get_terrain_weight(aStage, Vector2{(float)x, (float)y});

            if(nd + weight == d){
                d = nd, y = ny, x = nx;
                break;
            }
        }
    }
    path.push_back(Vector2{(float)(x+0.5), (float)(y+0.5)});
}

/// ウサギが次にジャンプすべき座標を計算
Vector2 next_path(const Stage& aStage, Path path, int& path_idx, Vector2 cur_pos, float power){
    Vector2 target;
    Terrain tar_terrain;
    Vector2 pre_point = cur_pos;
    Terrain pre_terrain = aStage.terrain(pre_point);
    bool first = true;
    int tmp_idx = 0;

    for(int k = std::min(path_idx + 5, (int)path.size()-1); k > path_idx; k--){
        auto tmp_target = path[k];
        if(same(tmp_target, cur_pos)){
            path_idx = k;
            break;
        }
    }

    while(1){
        // TODO : 逆戻り解消
        target = path[path_idx + tmp_idx];
        tar_terrain = aStage.terrain(target);

        auto next = aStage.getNextPos(cur_pos, power, target);
        auto nex_terrain = aStage.terrain(next);

        if(!first && tar_terrain > pre_terrain){
            target = pre_point;
            break;
        }
        //else if(path_idx+1 < (int)path.size() && same_float(target, next)){
        else if(path_idx+tmp_idx+1 < (int)path.size() && same_float(target, next)){
            path_idx++;
        }
        //else if(path_idx+1 < (int)path.size() && same(target, next)){
        else if(path_idx+tmp_idx+1 < (int)path.size() && same(target, next)){
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
}

/// TODO : {U, R, D, L}から出発した際の
/// 点間移動する際のウサギの通る座標 {rabbit, scroll01, scroll02, ...}
//Path sequence[4][MAX_V][MAX_V];

/// 点間の中心を移動する際のウサギの通る座標 {rabbit, scroll01, scroll02, ...}
Path sequence_simple[MAX_V][MAX_V];

/// 点間のダイクストラ距離 {rabbit, scroll01, scroll02, ...}
Path dijkstra_path[MAX_V][MAX_V];

/// dijkstra_pathを全点間で作成
void build_all_dijkstra_path(const Stage& aStage){
    const int scroll_num = aStage.scrolls().count();
    for(int i = 0; i < scroll_num+1; i++){
        for(int j = 1; j < scroll_num+1; j++){
            if(i == j) continue;
            auto src = (i == 0) ? aStage.rabbit().pos() : aStage.scrolls()[i-1];
            dijkstra_path[i][j].clear();
            //get_path_from_dijkstra(aStage, dijkstra_path[i][j], src.pos().y, src.pos().x, -1, j-1, -1);
            get_path_from_dijkstra_simple(aStage, dijkstra_path[i][j], src.pos().y, src.pos().x, -1, j-1, -1);
        }
    }
}

// v1 -> v2 へのウサギの通る座標(sequence_simple)を更新
void update_path_simple(const Stage& aStage, int v1, int v2, float power){
    auto scrolls = aStage.scrolls();
    auto dest = scrolls[v2-1].pos();
    auto pos = (v1 == 0) ? aStage.rabbit().pos() : scrolls[v1-1].pos();

    int path_idx = 0;
    sequence_simple[v1][v2].clear();
    while(1){
        auto target = next_path(aStage, dijkstra_path[v1][v2], path_idx, pos, power);
        pos = aStage.getNextPos(pos, power, target);
        sequence_simple[v1][v2].push_back(pos);
        //if(same(pos, dest))
        if(same_float(pos, dest))
            break;
    }
}

/// targets配列の初期化 & 
void targets_init(const Stage& aStage){
    const int nscrolls = aStage.scrolls().count();
    targets.resize(nscrolls);
    for(int i = 0; i < nscrolls; i++) 
        targets[i] = i;
}

void targets_shuffle(const Stage& aStage){
    const int nscrolls = aStage.scrolls().count();
    for(int i = 0; i < nscrolls; i++){
        int j = randxor() % nscrolls;
        std::swap(targets[i], targets[j]);
    }
}

/// targetsからpathをupdate
void update_path_from_targets(const Stage& aStage){
    const int nscrolls = aStage.scrolls().count();
    int src = 0;
    for(int i = 0; i < nscrolls; i++){
        auto idx = targets[i]+1;
        update_path_simple(aStage, src, idx, POWER[i]);
        src = idx;
    }
}

/// TODO : 部分的な高速計算？
/// 終了までに要するターン数を計算
int path_length(const Stage& aStage){
    const int nscrolls = aStage.scrolls().count();
    int src = 0;
    int npath = 0;
    for(int i = 0; i < nscrolls; i++){
        auto idx = targets[i] + 1;
        npath += sequence_simple[src][idx].size();
        src = idx;
    }
    return npath;
}

/// 終了までに要するターン数を計算
float path_length_from_dijkstra(const Stage& aStage){
    const auto scrolls = aStage.scrolls();
    const int nscrolls = scrolls.count();
    int src = -1;
    float npath = 0;
    for(int i = 0; i < nscrolls; i++){
        //auto idx = targets[i] + 1;
        auto dest = targets[i];
        //TODO : not size
        //npath += (float)dijkstra_path[src][idx].size() / POWER[i];
        Vector2 pos = {-1, -1};
        if(i == 0)
            pos = aStage.rabbit().pos();
        else
            pos = scrolls[src].pos();
        npath += (float)dist_scroll[dest][(int)pos.y][(int)pos.x] / (POWER[i]);
        src = dest;
    }
    return npath;
}

/// targetsのv1番目とv2番目を交換し、再計算
void change_vertex(const Stage& aStage, int v1, int v2){
    assert(v1 != v2);
    if(v1 > v2) std::swap(v1, v2);
    const int nscrolls = aStage.scrolls().count();
    int pre_v1 = (v1 == 0) ? 0 : targets[v1-1]+1;
    int pre_v2 = targets[v2-1]+1;
    std::swap(targets[v1], targets[v2]);
    if(v2 - v1 > 1){
        update_path_simple(aStage, pre_v1,        targets[v1]+1, POWER[v1]);
        update_path_simple(aStage, targets[v1]+1, targets[v1+1]+1, POWER[v1+1]);
        update_path_simple(aStage, pre_v2,        targets[v2]+1, POWER[v2]);
        if(v2+1 < nscrolls)
            update_path_simple(aStage, targets[v2]+1, targets[v2+1]+1, POWER[v2+1]);
    }
    else{
        update_path_simple(aStage, pre_v1,        targets[v1]+1, POWER[v1]);
        update_path_simple(aStage, targets[v1]+1, targets[v2]+1, POWER[v1+1]);
        if(v2+1 < nscrolls)
            update_path_simple(aStage, targets[v2]+1, targets[v2+1]+1, POWER[v2+1]);
    }
}

/// targetsのv1番目とv2番目を交換し、再計算
void change_vertex_simple(const Stage& aStage, int v1, int v2){
    assert(v1 != v2);
    std::swap(targets[v1], targets[v2]);
}

/// 2-opt法でTSPを解く (1.404sec)
void tsp_2opt(const Stage& aStage){
    const int nscrolls = aStage.scrolls().count();
    if(nscrolls > 1){
        int nstep = path_length(aStage);
        while(1){
            bool update = false;
            for(int i = 0; i < nscrolls; i++){
                for(int j = i+1; j < nscrolls; j++){
                    change_vertex(aStage, i, j);
                    int new_nstep = path_length(aStage);
                    if(new_nstep < nstep){
                        nstep = new_nstep;
                        update = true;
                    }
                    else{
                        //TODO : efficeint reverse
                        change_vertex(aStage, i, j);
                    }
                }
            }

            if(!update)
                break;
            //else
            //    printf("return %d\n", nstep);
        }
    }
}

/// 2-opt法でTSPを解く (1.404sec)
void tsp_2opt_simple(const Stage& aStage){
    const int nscrolls = aStage.scrolls().count();
    if(nscrolls > 1){
        float nstep = path_length_from_dijkstra(aStage);
        while(1){
            bool update = false;
            for(int i = 0; i < nscrolls; i++){
                for(int j = i+1; j < nscrolls; j++){
                    change_vertex_simple(aStage, i, j);
                    float new_nstep = path_length_from_dijkstra(aStage);
                    if(new_nstep < nstep){
                        nstep = new_nstep;
                        update = true;
                    }
                    else{
                        change_vertex_simple(aStage, i, j);
                    }
                }
            }

            if(!update)
                break;
            //else
            //    printf("return %d\n", nstep);
        }
    }
}

unsigned int randxor()
{
    static unsigned int x=123456789, y=362436069, z=521288629, w=88675123;
    unsigned int t;
    t = (x^(x<<11));
    x = y;
    y = z;
    z = w; 
    return( w=(w^(w>>19))^(t^(t>>8)) );
}

/// SA法でTSPを解く
void tsp_sa(const Stage& aStage, int iteration){
    const int N = aStage.scrolls().count();

    const double startTemp = 10;
    const double endTemp = 1;
    const int R = 10000;
    const int T = iteration;

    int nstep = path_length(aStage);
    for(int t = 0; t < iteration; t++){
        int v1 = randxor() % N, v2 = randxor() % N;
        while(v1 == v2)
            v1 = randxor() % N, v2 = randxor() % N;

        change_vertex(aStage, v1, v2);
        int new_nstep = path_length(aStage);

        double temp = startTemp + (endTemp - startTemp) * t / T;
        double probability = exp((nstep - new_nstep) / temp);
        bool force_next = probability > (double)(randxor() % R) / R;

        if(new_nstep < nstep || force_next){
            nstep = new_nstep;
        }
        else{
            //TODO : efficeint reverse
            change_vertex(aStage, v1, v2);
        }
    }
}

/// SA法でTSPを解く
void tsp_sa_simple(const Stage& aStage, int iteration){
    const int N = aStage.scrolls().count();

    const double startTemp = 100;
    const double endTemp = 1;
    const int R = 100000;
    const int T = iteration;

    float nstep = path_length(aStage);
    for(int t = 0; t < iteration; t++){
        int v1 = randxor() % N, v2 = randxor() % N;
        while(v1 == v2)
            v1 = randxor() % N, v2 = randxor() % N;

        change_vertex_simple(aStage, v1, v2);
        float new_nstep = path_length_from_dijkstra(aStage);

        double temp = startTemp + (endTemp - startTemp) * t / T;
        double probability = exp((nstep - new_nstep) / temp);
        bool force_next = probability > (double)(randxor() % R) / R;

        if(new_nstep < nstep || force_next){
            nstep = new_nstep;
        }
        else{
            change_vertex_simple(aStage, v1, v2);
        }
    }
}

Vector2 MygetTargetPos(const Stage& aStage){
    static int path_idx;
    static Path path;
    auto pos = aStage.rabbit().pos();
    if(aStage.turn() == 0 || same(pos, path[path.size()-1])){
        auto scrolls = aStage.scrolls();
        for(int i = 0; i < (int)targets.size(); i++){
            int idx = targets[i];
            auto scroll = scrolls[idx];
            if (!scroll.isGotten()) {
                path.clear();
                if(i == (int)targets.size()-1)
                    get_path_from_dijkstra(aStage, path, pos.y, pos.x, -1, idx, -1);
                else
                    get_path_from_dijkstra(aStage, path, pos.y, pos.x, -1, idx, -1);
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
    while(1){
        target = path[path_idx];
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
        else{
            if(nex_terrain == tar_terrain)
                target = next;
            break;
        }

        pre_point = target;
        pre_terrain = tar_terrain;
        first = false;
    }
    return target;
}


//------------------------------------------------------------------------------
/// コンストラクタ
/// @detail 最初のステージ開始前に実行したい処理があればここに書きます
Answer::Answer()
{
    float power = 1.0;
    for(int i = 0; i < 20; i++){
        POWER[i] = power;
        power *= 1.1;
    }
}

//------------------------------------------------------------------------------
/// デストラクタ
/// @detail 最後のステージ終了後に実行したい処理があればここに書きます
Answer::~Answer()
{
}

// for getTragetPos
//Path path;
//int path_idx = 0;
//bool big_jump;

/// TODO : Game Emulation
/// targets配列の順序をランダムに変えてみる？

//------------------------------------------------------------------------------
/// 各ステージ開始時に呼び出される処理
/// @detail 各ステージに対する初期化処理が必要ならここに書きます
/// @param aStage 現在のステージ
int test_idx;
void Answer::initialize(const Stage& aStage)
{
    test_idx++;
    //printf("test : %d\n", test_idx);

    // {rabbit, scroll0, scroll1, ...}
    Ver vertices;
    vertices_init(vertices, aStage);

    make_graph(aStage);
    for(int v = 0; v < (int)aStage.scrolls().count(); v++)
        dijkstra(aStage, v);
    calc_distance_dijkstra(aStage, vertices);

    build_all_dijkstra_path(aStage);

    //calc_distance_greedy(vertices);
    //const std::vector<double> Theta = {M_PI*sqrt(3.0/4), M_PI/6, 0};
    //const int topk = 10;
    //calc_distance_search_topk(aStage, vertices, topk, Theta);

    //sales_man_init(vertices.size());
    //sales_man_scroll(1 << 0, 0, vertices.size());
    //build_target_sequence_from_tsp();
    //update_path_from_targets(aStage);
    //tsp_2opt(aStage);

    const int nscrolls = aStage.scrolls().count();
    targets_init(aStage);
    if(nscrolls >= 3){
        //sales_man_init(vertices.size());
        //sales_man_scroll(1 << 0, 0, vertices.size());
        //build_target_sequence_from_tsp();
        //update_path_from_targets(aStage);
    
        //const int iteration = 2500000;
        //const int iteration = 50000;
        // 十分っぽい
        const int iteration = 50000;
        std::vector<int> best_targets;
        int best_score = 100000;
        //const int nloop = (nscrolls < 8) ? 10 : 30;
        const int nloop = (nscrolls < 8) ? 10 : 30;
        for(int i = 0; i < nloop; i++){
            targets_shuffle(aStage);

            tsp_sa_simple(aStage, iteration);
            tsp_2opt_simple(aStage);

            float dist = path_length_from_dijkstra(aStage);
            std::vector<int> _targets = targets;
            tsp_2opt(aStage);
            float new_dist = path_length_from_dijkstra(aStage);
            if(new_dist > dist)
                targets = _targets;

            //TODO : test
            EmulateGame eGame(aStage);
            int score = eGame.run();
            if(score < best_score){
                best_score = score;
                best_targets = targets;
            }
        }

        // TODO : efficient targets copy
        targets = best_targets;

    }
    else{
        float min_dist = 1e10;
        std::vector<int> _targets;
        do{
            float new_dist = path_length_from_dijkstra(aStage);
            if(new_dist < min_dist){
                _targets = targets;
                min_dist = new_dist;
            }
        }while(std::next_permutation(targets.begin(), targets.end()));
        targets = _targets;
    }
}

//------------------------------------------------------------------------------
/// 毎フレーム呼び出される処理
/// @detail 移動先を決定して返します
/// @param aStage 現在のステージ
/// @return 移動の目標座標
Vector2 Answer::getTargetPos(const Stage& aStage)
{
    return MygetTargetPos(aStage);
/*
    auto pos = aStage.rabbit().pos();
    for(auto scroll : aStage.scrolls()) {
        // まだ手に入れていない巻物を探して、そこに向かって飛ぶ
        if (!scroll.isGotten()) {
            return scroll.pos();
        }
    }
    return pos;
*/
/*
    auto pos = aStage.rabbit().pos();
    static int path_idx;
    if(aStage.turn() == 0 || same(pos, path[path.size()-1])){
        auto scrolls = aStage.scrolls();
        for(int i = 0; i < (int)targets.size(); i++){
            int idx = targets[i];
            auto scroll = scrolls[idx];
            if (!scroll.isGotten()) {
                int src = (i == 0) ? -1 : i-1;
                path.clear();
                if(i == (int)targets.size()-1)
                    get_path_from_dijkstra(aStage, path, pos.y, pos.x, src, idx, -1);
                else
                    //get_path_from_dijkstra(aStage, path, pos.y, pos.x, src, idx, targets[i+1]);
                    get_path_from_dijkstra(aStage, path, pos.y, pos.x, src, idx, -1);
                path_idx = 0;
                break;
            }
        }
    }
    auto next = next_path(aStage, path, path_idx, pos, aStage.rabbit().power());
    return next;
*/
    //TODO : 池を後の方に回す
/*
    static int src, dest;
    static int dest_idx;
    static int path_idx;
    if(aStage.turn() == 0){
        dest_idx = 0;
        src = 0;
        dest = 1 + targets[dest_idx++];
        path_idx = 0;
    }

    if((int)sequence_simple[src][dest].size() == path_idx){
        //printf("Change target\n\n");
        src = dest;
        dest = 1 + targets[dest_idx++];
        path_idx = 0;
    }

    auto next = sequence_simple[src][dest][path_idx++];
    return next;
*/
/*
    static int path_idx;
    static Path path;
    auto pos = aStage.rabbit().pos();
    if(aStage.turn() == 0 || same(pos, path[path.size()-1])){
        auto scrolls = aStage.scrolls();
        for(int i = 0; i < (int)targets.size(); i++){
            int idx = targets[i];
            auto scroll = scrolls[idx];
            if (!scroll.isGotten()) {
                path.clear();
                if(i == (int)targets.size()-1)
                    get_path_from_dijkstra(aStage, path, pos.y, pos.x, -1, idx, -1);
                else
                    get_path_from_dijkstra(aStage, path, pos.y, pos.x, -1, idx, -1);
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
    while(1){
        target = path[path_idx];
        tar_terrain = aStage.terrain(target);
 
        auto next = aStage.getNextPos(pos, aStage.rabbit().power(), target);
        auto nex_terrain = aStage.terrain(next);
 
        if(!first && tar_terrain > pre_terrain){
            /// TODO:pre_pointがrabbit.pos()の場合を除く。
            target = pre_point;
            break;
        }
        else if(path_idx+1 < (int)path.size() && same_float(target, next)){
            path_idx++;
        }
        else{
            if(nex_terrain == tar_terrain)
                target = next;
            break;
        }
 
        pre_point = target;
        pre_terrain = tar_terrain;
        first = false;
    }
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
