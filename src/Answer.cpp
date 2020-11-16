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
typedef std::pair<int, int> P;

const float EPS = 1e-5;
const double INF = 1e50;
/// rabbit : 0, others : 1 - 20
const int MAX_V = 1 + Parameter::MaxScrollCount;

/// 巻物をn個とった際のジャンプ力
float POWER[30];

/// ウサギ&巻物 間の速さをジャンプ力を1.0とした場合の移動時間
double distance[25][25];

/// dijkstra方で用いる距離
int dist_scroll[Parameter::MaxScrollCount][Parameter::StageHeight][Parameter::StageWidth];

/// 巻物を取る順番
std::vector<int> targets;

struct Edge{int y, x, weight;};

/// 各タイルを結んだグラフ {nx, ny, weight}
std::vector<Edge> Graph[Parameter::StageHeight][Parameter::StageWidth];

void calc_distance_dijkstra(const Stage& aStage, Ver vertices);
void sales_man_init(int size);
double sales_man(int S, int v, int size);
void vertices_init(Ver& vertices, const Stage& aStage);
void build_target_sequence_from_tsp();
void make_graph(const Stage& aStage);
bool same_float(Vector2 v1, Vector2 v2);
bool same(Vector2 v1, Vector2 v2);
int get_terrain_weight(const Stage& aStage, Vector2 pos);

unsigned int randxor();
void get_path_from_dijkstra(const Stage& aStage, Path& path, int y1, int x1, int src, int dest, int dest2, bool insert_first, bool init_dir, bool set_dest_center);
/// {y1, x1} から scrolls[dest] までの経路を復元. 実装を簡単にするため、目標値はすべてマスの中心にする
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
Vector2 MygetTargetPos(const Stage& aStage, int n);


////////////////////////////////////////////////////////////////////////

class EmulateGame{
    private:
        Stage stage;

    public:
        EmulateGame(const Stage& aStage){
            init(aStage);
        }

        void init(const Stage& aStage){
            stage = aStage;
        }

        int run(int n){
            int current_turn = 0;
            while (!stage.isEnd() && stage.turn() < Parameter::GameTurnLimit && !stage.isOutOfBounds(stage.rabbit().pos())) {
                // ターン開始
                auto targetPos = MygetTargetPos(stage, n);
                stage.update(targetPos);
                stage.advanceTurn();
                current_turn++;
            }
            return current_turn;
        }
};

/// pos地点の地形から設定した重みを返す
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

/// v1とv2が同じ位置かどうか判定
bool same_float(Vector2 v1, Vector2 v2){
    return (fabs(v1.x - v2.x) < EPS && fabs(v1.y - v2.y) < EPS);
}

/// v1とv2が同じマスかどうか判定
bool same(Vector2 v1, Vector2 v2){
    //return ((int)v1.x == (int)v2.x && (int)v1.y == (int)v2.y);
    return (static_cast<int>(v1.x) == static_cast<int>(v2.x) 
            && static_cast<int>(v1.y) == static_cast<int>(v2.y));
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
/// TODO : 何らかのパラメータでアンサンブル
/// 最後の地点の設定
/// 4隅をうまく使えないか？
/// TODO : 最後の地点ではなく、次に曲がる方向に寄っておく
Vector2 first_step[MAX_V];
void get_path_from_dijkstra(const Stage& aStage, Path& path, int y1, int x1, int src, int dest, int dest2, bool insert_first, bool init_dir, bool set_dest_center){
    assert(src == -1);
    //TODO : 最初の地点を保存？
    int d = dist_scroll[dest][y1][x1];
    auto target = aStage.scrolls()[dest].pos();
    int y = y1, x = x1;
    if(insert_first)
        path.push_back(Vector2{(float)(x+0.5), (float)(y+0.5)});
    while(d != 0){
        //TODO : to center? 
        // 特にtargetが池の場合は少し入ってすぐ出るべきだから中心はまずい
        //path.push_back(Vector2{(float)x+(float)0.5, (float)y+(float)0.5});
        int dx[2][4][4] = {{{1, 0, -1, 0}, {0, -1, 0, 1}, {-1, 0, 1, 0}, {0, 1, 0, -1}},
                           {{0, -1, 0, 1}, {-1, 0, 1, 0}, {0, 1, 0, -1}, {1, 0, -1, 0}}};
        int dy[2][4][4] = {{{0, 1, 0, -1}, {1, 0, -1, 0}, {0, -1, 0, 1}, {-1, 0, 1, 0}},
                           {{1, 0, -1, 0}, {0, -1, 0, 1}, {-1, 0, 1, 0}, {0, 1, 0, -1}}};
        int direction = get_direction(Vector2{(float)(x+0.5), (float)(y+0.5)}, target);
        const int turn = path.size() + init_dir;

        for(int i = 0; i < 4; i++){
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
                        path.push_back(Vector2{tx, ty});
                        d = nd, y = ny, x = nx;
                        break;
                    }
                    else{
                        tx = ((1+EPS)*(nx+0.5) + (x+0.5)) / (2. + EPS);
                        ty = ((1+EPS)*(ny+0.5) + (y+0.5)) / (2. + EPS);
                        path.push_back(Vector2{tx, ty});

                        if(set_dest_center){
                            tx = (x + nx + 1) / 2.0;
                            ty = (y + ny + 1) / 2.0;
                            path.push_back(Vector2{tx, ty});
                        }

                        //// dest から dest2の方向の隅を最終到着点に設定
                        auto nex_target = aStage.scrolls()[dest2].pos();
                        int dir = get_direction(target, nex_target);
                        // upper-right
                        if(dir == 0){
                            tx = nx + 1 - EPS; 
                            ty = ny + 1 - EPS; 
                        }
                        // upper-left
                        else if(dir == 1){
                            tx = nx; 
                            ty = ny + 1 - EPS; 
                        }
                        // lower-left
                        else if(dir == 2){
                            tx = nx;
                            ty = ny;
                        }
                        // lower-right
                        else{
                            tx = nx + 1 - EPS;
                            ty = ny;
                        }
                        path.push_back(Vector2{tx, ty});
                        d = nd, y = ny, x = nx;
                        break;
                    }
                }else{
                    auto ter_org = aStage.terrain(Vector2{x+EPS, y+EPS});
                    auto ter_nxt = aStage.terrain(Vector2{nx+EPS, ny+EPS});
                    if(ter_org == ter_nxt){
                        tx = (x + nx + 1) / 2.0;
                        ty = (y + ny + 1) / 2.0;
                        //if(x == nx){
                        //    ty = (y + ny + 1) / 2.0;
                        //    if((int)target.x == x)
                        //        tx = (x + nx + 1) / 2.0;
                        //    else if((int)target.x > x)
                        //        tx = x + 1 - EPS;
                        //    else
                        //        tx = x;
                        //}
                        //else{
                        //    tx = (x + nx + 1) / 2.0;
                        //    if((int)target.y == y)
                        //        ty = (y + ny + 1) / 2.0;
                        //    else if((int)target.y > y)
                        //        ty = y + 1 - EPS;
                        //    else
                        //        ty = y;
                        //}
                    }
                    else if(ter_org < ter_nxt){
                        //// TODO : ensemble???
                        tx = ((nx+0.5) + (1+EPS)*(x+0.5)) / (2. + EPS);
                        ty = ((ny+0.5) + (1+EPS)*(y+0.5)) / (2. + EPS);
                        //if(x == nx){
                        //    ty = ((ny+0.5) + (1+EPS)*(y+0.5)) / (2. + EPS);
                        //    if((int)target.x == x)
                        //        tx = (x + nx + 1) / 2.0;
                        //    else if((int)target.x > x)
                        //        tx = x + 1 - EPS;
                        //    else
                        //        tx = x;
                        //}
                        //else{
                        //    tx = ((nx+0.5) + (1+EPS)*(x+0.5)) / (2. + EPS);
                        //    if((int)target.y == y)
                        //        ty = (y + ny + 1) / 2.0;
                        //    else if((int)target.y > y)
                        //        ty = y + 1 - EPS;
                        //    else
                        //        ty = y;
                        //}
                    }
                    else{
                        tx = ((1+EPS)*(nx+0.5) + (x+0.5)) / (2. + EPS);
                        ty = ((1+EPS)*(ny+0.5) + (y+0.5)) / (2. + EPS);
                        //if(x == nx){
                        //    ty = ((1+EPS)*(ny+0.5) + (y+0.5)) / (2. + EPS);
                        //    if((int)target.x == x)
                        //        tx = (x + nx + 1) / 2.0;
                        //    else if((int)target.x > x)
                        //        tx = x + 1 - EPS;
                        //    else
                        //        tx = x;
                        //}
                        //else{
                        //    tx = ((1+EPS)*(nx+0.5) + (x+0.5)) / (2. + EPS);
                        //    if((int)target.y == y)
                        //        ty = (y + ny + 1) / 2.0;
                        //    else if((int)target.y > y)
                        //        ty = y + 1 - EPS;
                        //    else
                        //        ty = y;
                        //}
                    }

                    if(path.size() == 0 && src >= 0)
                        first_step[src] = Vector2{tx, ty};

                    path.push_back(Vector2{tx, ty});
                    d = nd, y = ny, x = nx;
                    break;
                }
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

/// xor128
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

Vector2 Answer01(const Stage& aStage, bool path_init_dir, bool path_dest_center, bool strict){
    static int path_idx;
    static Path path;
    auto pos = aStage.rabbit().pos();
    /// targets配列の順番毎にpathに移動すべき点を格納
    if(aStage.turn() == 0 || same(pos, path[path.size()-1])){
        auto scrolls = aStage.scrolls();
        for(int i = 0; i < (int)targets.size(); i++){
            int idx = targets[i];
            auto scroll = scrolls[idx];
            if (!scroll.isGotten()) {
                path.clear();
                int nex_dest = (i+1 < (int)targets.size()) ? targets[i+1] : -1;
                if(i == (int)targets.size()-1)
                    //get_path_from_dijkstra(aStage, path, pos.y, pos.x, -1, idx, -1, false, path_init_dir);
                    get_path_from_dijkstra(aStage, path, pos.y, pos.x, -1, idx, nex_dest, false, path_init_dir, path_dest_center);
                else
                    //get_path_from_dijkstra(aStage, path, pos.y, pos.x, -1, idx, -1, false, path_init_dir);
                    get_path_from_dijkstra(aStage, path, pos.y, pos.x, -1, idx, nex_dest, false, path_init_dir, path_dest_center);
                path_idx = 0;
                break;
            }
        }
    }

/// 到達可能な点の次の点を目標とする。
    Vector2 target;
    Terrain tar_terrain;
    Vector2 pre_point = pos;
    Terrain pre_terrain = aStage.terrain(pre_point);
    bool first = true;
    while(1){
        target = path[path_idx];
        tar_terrain = aStage.terrain(target);

        auto next = aStage.getNextPos(pos, aStage.rabbit().power(), target);
        //auto nex_terrain = aStage.terrain(next);

        if(!first && tar_terrain > pre_terrain){
            target = pre_point;
            break;
        }
        else if(path_idx+1 < (int)path.size() && ((strict && same_float(target, next)) || (!strict && same(target, next)))){
            path_idx++;
        }
        else{
            break;
        }

        pre_point = target;
        pre_terrain = tar_terrain;
        first = false;
    }
    return target;
}

Vector2 Answer04(const Stage& aStage, int postk, bool path_init_dir, bool path_dest_center){
    static int path_idx;
    static Path path;
    auto pos = aStage.rabbit().pos();
    static std::map<std::pair<int, int>, int> coord2idx;
    static bool jump;
    /// targets配列の順番毎にpathに移動すべき点を格納
    if(aStage.turn() == 0 || same(pos, path[path.size()-1])){
        auto scrolls = aStage.scrolls();
        for(int i = 0; i < (int)targets.size(); i++){
            int idx = targets[i];
            auto scroll = scrolls[idx];
            if (!scroll.isGotten()) {
                path.clear();
                int nex_dest = (i+1 < (int)targets.size()) ? targets[i+1] : -1;
                if(i == (int)targets.size()-1)
                    get_path_from_dijkstra(aStage, path, pos.y, pos.x, -1, idx, nex_dest, true, path_init_dir, path_dest_center);
                else
                    get_path_from_dijkstra(aStage, path, pos.y, pos.x, -1, idx, nex_dest, true, path_init_dir, path_dest_center);
                path_idx = 0;
                jump = false;
                coord2idx.clear();
                for(int k = 0; k < (int)path.size(); k++){
                    auto coord = path[k];
                    coord2idx[std::make_pair((int)coord.x, (int)coord.y)] = k;
                }
                break;
            }
        }
    }

    auto p = std::make_pair((int)pos.x, (int)pos.y);
    static Vector2 checkpoint;
    /// postk個先の点を目標とする。その際地形が悪化する場合は前の点を目標とする。
    if(coord2idx.find(p) == coord2idx.end()){
        return checkpoint;
    }
    else{
        int cur_path_idx = coord2idx[p];
        auto cur_terrain = aStage.terrain(pos);
        // TODO merge
        for(int k = postk; k >= 0; k--){
            if(cur_path_idx + k >= (int)path.size())
                continue;
            int next_path_idx = cur_path_idx + k;
            auto target = path[next_path_idx];
            auto tmp_pos = pos;
            bool valid = true;
            while(1){
                auto next = aStage.getNextPos(tmp_pos, aStage.rabbit().power(), target);
                auto nex_terrain = aStage.terrain(next);
                if(jump){
                    jump = false;
                    cur_terrain = nex_terrain;
                }
                if(same(next, target)){
                    break;
                }
                else if(nex_terrain > cur_terrain){
                    valid = false;
                    break;
                }
                tmp_pos = next;
            }

            if(valid){
                path_idx = next_path_idx;
                jump = (k == 0);
                break;
            }
        }

        checkpoint = path[path_idx];
        return path[path_idx];
    }

}

/// 次の地点
/// TODO : 地点毎にアンサンブル
//#define NANSWER 14
#define NANSWER 28
//#define NANSWER 1
Vector2 MygetTargetPos(const Stage& aStage, int n){
    assert(0 <= n && n < NANSWER);
    switch(n){
        case 0 : return Answer01(aStage, true, false, true);
        case 1 : return Answer01(aStage, false, false, true);
        case 2 : return Answer04(aStage, 3, false, false);
        case 3 : return Answer04(aStage, 4, false, false);
        case 4 : return Answer04(aStage, 5, false, false);
        case 5 : return Answer04(aStage, 8, false, false);
        case 6 : return Answer04(aStage, 12, false, false);
        case 7 : return Answer04(aStage, 3, true, false);
        case 8 : return Answer04(aStage, 4, true, false);
        case 9 : return Answer04(aStage, 5, true, false);
        case 10 : return Answer04(aStage, 8, true, false);
        case 11 : return Answer04(aStage, 12, true, false);
        case 12 : return Answer01(aStage, true, false, false);
        case 13 : return Answer01(aStage, false, false, false);

        case 14 : return Answer01(aStage, true, true, true);
        case 15 : return Answer01(aStage, false, true, true);
        case 16 : return Answer04(aStage, 3, false, true);
        case 17 : return Answer04(aStage, 4, false, true);
        case 18 : return Answer04(aStage, 5, false, true);
        case 19 : return Answer04(aStage, 8, false, true);
        case 20 : return Answer04(aStage, 12, false, true);
        case 21 : return Answer04(aStage, 3, true, true);
        case 22 : return Answer04(aStage, 4, true, true);
        case 23 : return Answer04(aStage, 5, true, true);
        case 24 : return Answer04(aStage, 8, true, true);
        case 25 : return Answer04(aStage, 12, true, true);
        case 26 : return Answer01(aStage, true, true, false);
        case 27 : return Answer01(aStage, false, true, false);
    }
    assert(true);
    return Answer01(aStage, true, true, true);
}

Vector2 execute_answer(const Stage& aStage, int nth_answer, Path path, Vector2 cur_pos, float cur_power, bool init, bool strict, int postk){
    static int path_idx;
    if(init) {
        path_idx = 0;
    }

    //auto pos = aStage.rabbit().pos();
    auto pos = cur_pos;
    if(nth_answer == 1){
        /// 到達可能な点の次の点を目標とする。
        Vector2 target;
        Terrain tar_terrain;
        Vector2 pre_point = pos;
        Terrain pre_terrain = aStage.terrain(pre_point);
        bool first = true;
        while(1){
            target = path[path_idx];
            tar_terrain = aStage.terrain(target);

            auto next = aStage.getNextPos(pos, cur_power, target);
            //auto nex_terrain = aStage.terrain(next);

            if(!first && tar_terrain > pre_terrain){
                target = pre_point;
                break;
            }
            else if(path_idx+1 < (int)path.size() && ((strict && same_float(target, next)) || (!strict && same(target, next)))){
                path_idx++;
            }
            else{
                break;
            }

            pre_point = target;
            pre_terrain = tar_terrain;
            first = false;
        }
        return target;
    }
    else if(nth_answer == 4){
        static std::map<std::pair<int, int>, int> coord2idx;
        static int jump;
        if(init){
            jump = 0;
            coord2idx.clear();
            for(int k = 0; k < (int)path.size(); k++){
                auto coord = path[k];
                coord2idx[std::make_pair((int)coord.x, (int)coord.y)] = k;
            }
        }

        auto p = std::make_pair((int)pos.x, (int)pos.y);
        static Vector2 checkpoint;
        /// postk個先の点を目標とする。その際地形が悪化する場合は前の点を目標とする。
        if(coord2idx.find(p) == coord2idx.end()){
            return checkpoint;
        }
        else{
            int cur_path_idx = coord2idx[p];
            auto cur_terrain = aStage.terrain(pos);
            // TODO merge
            for(int k = postk; k >= 0; k--){
                if(cur_path_idx + k >= (int)path.size())
                    continue;
                int next_path_idx = cur_path_idx + k;
                auto target = path[next_path_idx];
                auto tmp_pos = pos;
                bool valid = true;
                while(1){
                    auto next = aStage.getNextPos(tmp_pos, cur_power, target);
                    auto nex_terrain = aStage.terrain(next);
                    if(jump){
                        jump = false;
                        cur_terrain = nex_terrain;
                    }
                    if(same(next, target)){
                        break;
                    }
                    else if(nex_terrain > cur_terrain){
                        valid = false;
                        break;
                    }
                    tmp_pos = next;
                }

                if(valid){
                    path_idx = next_path_idx;
                    jump = (k == 0);
                    break;
                }
            }

            checkpoint = path[path_idx];
            return path[path_idx];
        }
    }
    else{
        assert(false);
    }
    return Vector2{-1, -1};
}

void set_ensemble_parameters(int n, int& nth_answer, bool& strict, int& postk, bool& path_init_dir, bool& path_dest_center){
    assert(0 <= n && n < NANSWER);
    switch(n){
        case 0 :  nth_answer = 1, path_init_dir = true,  path_dest_center = false, strict = true; break;
        case 1 :  nth_answer = 1, path_init_dir = false, path_dest_center = false, strict = true; break;
        case 2 :  nth_answer = 4, postk = 3,  path_init_dir = false, path_dest_center = false; break;
        case 3 :  nth_answer = 4, postk = 4,  path_init_dir = false, path_dest_center = false; break;
        case 4 :  nth_answer = 4, postk = 5,  path_init_dir = false, path_dest_center = false; break;
        case 5 :  nth_answer = 4, postk = 8,  path_init_dir = false, path_dest_center = false; break;
        case 6 :  nth_answer = 4, postk = 12, path_init_dir = false, path_dest_center = false; break;
        case 7 :  nth_answer = 4, postk = 3,  path_init_dir = true,  path_dest_center = false; break;
        case 8 :  nth_answer = 4, postk = 4,  path_init_dir = true,  path_dest_center = false; break;
        case 9 :  nth_answer = 4, postk = 5,  path_init_dir = true,  path_dest_center = false; break;
        case 10 : nth_answer = 4, postk = 8,  path_init_dir = true,  path_dest_center = false; break;
        case 11 : nth_answer = 4, postk = 12, path_init_dir = true,  path_dest_center = false; break;
        case 12 : nth_answer = 1, path_init_dir = true,  path_dest_center = false, strict = false; break;
        case 13 : nth_answer = 1, path_init_dir = false, path_dest_center = false, strict = false; break;

        case 14 : nth_answer = 1, path_init_dir = true,  path_dest_center = true, strict = true; break;
        case 15 : nth_answer = 1, path_init_dir = false, path_dest_center = true, strict = true; break;
        case 16 : nth_answer = 4, postk = 3,  path_init_dir = false, path_dest_center = true; break;
        case 17 : nth_answer = 4, postk = 4,  path_init_dir = false, path_dest_center = true; break;
        case 18 : nth_answer = 4, postk = 5,  path_init_dir = false, path_dest_center = true; break;
        case 19 : nth_answer = 4, postk = 8,  path_init_dir = false, path_dest_center = true; break;
        case 20 : nth_answer = 4, postk = 12, path_init_dir = false, path_dest_center = true; break;
        case 21 : nth_answer = 4, postk = 3,  path_init_dir = true,  path_dest_center = true; break;
        case 22 : nth_answer = 4, postk = 4,  path_init_dir = true,  path_dest_center = true; break;
        case 23 : nth_answer = 4, postk = 5,  path_init_dir = true,  path_dest_center = true; break;
        case 24 : nth_answer = 4, postk = 8,  path_init_dir = true,  path_dest_center = true; break;
        case 25 : nth_answer = 4, postk = 12, path_init_dir = true,  path_dest_center = true; break;
        case 26 : nth_answer = 1, path_init_dir = true,  path_dest_center = true, strict = false; break;
        case 27 : nth_answer = 1, path_init_dir = false, path_dest_center = true, strict = false; break;
    }
}

int answers[30];
Path ans_path[50];
int ans_num;

/// 各answerで最もよいものをanswersに格納
void search_best_answers(const Stage& aStage, int nth_loop){
    ans_path[nth_loop].clear();
    auto scrolls = aStage.scrolls();
    auto cur_pos = aStage.rabbit().pos();
    //// 各地点間で最も良い方法を探索
    for(int i = 0; i < (int)targets.size(); i++){
        int dest = targets[i];
        int nex_dest = (i+1 < (int)targets.size()) ? targets[i+1] : -1;
        auto scroll = scrolls[dest];
        //auto cur_pos = (i == 0) ? aStage.rabbit().pos() : scrolls[targets[i-1]].pos();
        auto power = POWER[i];
        int best_turn = 10000;
        int best_answer = -1;
        Vector2 best_last_pos;

        /// 各方法での通ったパスを記録
        /// TODO : 計算時間
        Path checkpoints[NANSWER];

        //// 各answerを試す
        for(int n = 0; n < NANSWER; n++){
            auto tmp_pos = cur_pos;
            bool init = true;

            int nth_answer = 0, postk = 0;
            bool strict = false;
            bool path_init_dir = false; 
            bool path_dest_center = false;

            set_ensemble_parameters(n, nth_answer, strict, postk, path_init_dir, path_dest_center);

            ////TODO : skip if using same path
            Path path;
            get_path_from_dijkstra(aStage, path, tmp_pos.y, tmp_pos.x, -1, dest, nex_dest, nth_answer == 4, path_init_dir, path_dest_center);

            int current_turn = 0;
            //while (!stage.isEnd() && stage.turn() < Parameter::GameTurnLimit && !stage.isOutOfBounds(stage.rabbit().pos())) {
            while (!same(scroll.pos(), tmp_pos) && current_turn < Parameter::GameTurnLimit) {
                // ターン開始
                //auto targetPos = MygetTargetPos(stage, n);
                auto targetPos = execute_answer(aStage, nth_answer, path, tmp_pos, power, init, strict, postk);
                tmp_pos = aStage.getNextPos(tmp_pos, power, targetPos);
                checkpoints[n].push_back(tmp_pos);
                current_turn++;
                init = false;
            }
            //printf("Score %d : %d\n", n, current_turn);

            if(current_turn < best_turn){
                best_turn = current_turn;
                best_answer = n;
                best_last_pos = tmp_pos;
            }
        }
        //printf("best : %d\n", best_answer);
        answers[i] = best_answer;
        cur_pos = best_last_pos;
        for(auto&& v : checkpoints[best_answer])
            ans_path[nth_loop].push_back(v);
        //printf("path size : %d\n", (int)ans_path[nth_loop].size());
    }
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
/// for choosing path algo
int answer;

//------------------------------------------------------------------------------
/// 各ステージ開始時に呼び出される処理
/// @detail 各ステージに対する初期化処理が必要ならここに書きます
/// @param aStage 現在のステージ
void Answer::initialize(const Stage& aStage)
{
    //static int test_idx++;
    //printf("test : %d\n", test_idx);

    /// {rabbit, scroll0, scroll1, ...}のようにvertices配列にそれぞれの座標を格納
    Ver vertices;
    vertices_init(vertices, aStage);

    /// 各マス毎に{移動先、辺の重み}をGraph配列に格納
    make_graph(aStage);
    /// 各巻物からすべてのマスへの最短距離を求める
    for(int v = 0; v < (int)aStage.scrolls().count(); v++)
        dijkstra(aStage, v);
    /// ウサギ&巻物 間の距離をdistance配列に格納
    calc_distance_dijkstra(aStage, vertices);

    build_all_dijkstra_path(aStage);

    const int nscrolls = aStage.scrolls().count();
    targets_init(aStage);
    if(nscrolls >= 3){
        /*
        std::vector<int> best_targets;
        int best_score = 100000;
        //const int iteration = 2500000;
        // 十分っぽい
        //// 27307
        //const int iteration = 500000;
        //// 27340
        const int iteration = 50000;
        //const int nloop = (nscrolls < 8) ? 1 : 1;
        //// 28060
        const int nloop = (nscrolls < 8) ? 10 : 30;
        //// 28012
        //const int nloop = (nscrolls < 8) ? 30 : 90;
        /// nloop個回探索
        for(int i = 0; i < nloop; i++){
            targets_shuffle(aStage);

            /// SA法をなどを用いて探索
            tsp_sa_simple(aStage, iteration);
            tsp_2opt_simple(aStage);

            //// 1.8sec for 20 elms
            float dist = path_length_from_dijkstra(aStage);
            std::vector<int> _targets = targets;
            tsp_2opt(aStage);
            float new_dist = path_length_from_dijkstra(aStage);
            if(new_dist > dist)
                targets = _targets;

            /// NANSWER個の移動アルゴリズムから最も良いものを選択
            for(int k = 0; k < NANSWER; k++){
                EmulateGame eGame(aStage);
                int score = eGame.run(k);
                if(score < best_score){
                    best_score = score;
                    best_targets = targets;
                    answer = k;
                }
            }
        }
        targets = best_targets;
        */

        int best_score = 100000;
        const int iteration = 50000;
        const int nloop = (nscrolls < 8) ? 10 : 30;
        //const int nloop = (nscrolls < 8) ? 1 : 1;
        /// nloop個回探索
        for(int i = 0; i < nloop; i++){
            targets_shuffle(aStage);

            /// SA法をなどを用いて探索
            tsp_sa_simple(aStage, iteration);
            tsp_2opt_simple(aStage);

            //// 1.8sec for 20 elms
            float dist = path_length_from_dijkstra(aStage);
            std::vector<int> _targets = targets;
            tsp_2opt(aStage);
            float new_dist = path_length_from_dijkstra(aStage);
            if(new_dist > dist)
                targets = _targets;

            /// NANSWER個の移動アルゴリズムから最も良いものを選択
            search_best_answers(aStage, i);
            int score = ans_path[i].size();
            if(score < best_score){
                best_score = score;
                ans_num = i;
            }
        }
    }
    else{
        /*
        int best_score = 100000;
        std::vector<int> best_targets;
        do{
            for(int k = 0; k < NANSWER; k++){
                EmulateGame eGame(aStage);
                int score = eGame.run(k);
                if(score < best_score){
                    best_score = score;
                    best_targets = targets;
                    answer = k;
                }
            }
        }while(std::next_permutation(targets.begin(), targets.end()));
        targets = best_targets;
        */
        int best_score = 100000;
        int idx = 0;
        do{
            search_best_answers(aStage, idx);
            int score = ans_path[idx].size();
            if(score < best_score){
                best_score = score;
                ans_num = idx;
            }
            idx++;
        }while(std::next_permutation(targets.begin(), targets.end()));
    }
}

//------------------------------------------------------------------------------
/// 毎フレーム呼び出される処理
/// @detail 移動先を決定して返します
/// @param aStage 現在のステージ
/// @return 移動の目標座標
Vector2 Answer::getTargetPos(const Stage& aStage)
{
    //return MygetTargetPos(aStage, answer);

    static int path_idx;
    if(aStage.turn() == 0)
        path_idx = 0;
    return ans_path[ans_num][path_idx++];
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
