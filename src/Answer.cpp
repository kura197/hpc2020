//------------------------------------------------------------------------------
/// @file
/// @author   ハル研究所プログラミングコンテスト実行委員会
///
/// @copyright  (C)HAL Laboratory, Inc.
/// @attention  このファイルの利用は、同梱のREADMEにある
///             利用条件に従ってください。
//------------------------------------------------------------------------------

#include "Answer.hpp"
#include <vector>


//------------------------------------------------------------------------------
namespace hpc {

typedef std::vector<Vector2> Ver;

const double INF = 1e50;
/// rabbit : 0, others : 1 - 20
const int MAX_V = 21;

/// ウサギ&巻物 間の速さをジャンプ力を1.0とした場合の移動時間
double distance[25][25];

/// 巻物の取る順番
std::vector<int> targets;

double dist(Vector2 v1, Vector2 v2){
    double diff_y = v1.y - v2.y;
    double diff_x = v1.x - v2.x;
    return std::sqrt(diff_y*diff_y + diff_x*diff_x);
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
void target_sequence(){
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

//------------------------------------------------------------------------------
/// 各ステージ開始時に呼び出される処理
/// @detail 各ステージに対する初期化処理が必要ならここに書きます
/// @param aStage 現在のステージ
void Answer::initialize(const Stage& aStage)
{
    // {rabbit, scroll0, scroll1, ...}
    Ver vertices;
    vertices_init(vertices, aStage);

    calc_distance_greedy(vertices);

    sales_man_init(vertices.size());
    sales_man(1 << 0, 0, vertices.size());

    target_sequence();
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
*/
    auto scrolls = aStage.scrolls();
    for(int idx : targets){
        auto scroll = scrolls[idx];
        if (!scroll.isGotten()) {
            return scroll.pos();
        }
    }

    return pos;
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
