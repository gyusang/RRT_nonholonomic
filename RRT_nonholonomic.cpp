/* 
*    Author: Lee Sanggyu, SNU CSE 2020
*    email: effer@snu.ac.kr
*/

#include <iostream>
#include <cmath>
#include <ctime>
#include <vector>
#include <fstream>
#include <algorithm>

#define N 2000
#define PIXEL_WIDTH 10
#define PIXEL_HEIGHT 10
#define DENSITY 10 // pixels per meter
#define TERMINATE_DIST 0.05
#define NEAR_RADIUS 0.2
#define DT 0.02
#define MAX_V 1
#define MAX_OMEGA 5
#define GAMMA 0.05
#define MAX_STEER_ITER 1000
#define SAMPLE_GOAL_RATE 0.05
#define ANGLE_DISTANCE_FACTOR 0.1

using namespace std;

double randDouble(double min, double max){
    double randvar = (double)rand() / (double)RAND_MAX;
    return (max-min)*randvar + min;
}

class Obstacle {
    public:
        bool* occupancy[PIXEL_HEIGHT];
        bool checkFree(int u, int v);
        bool checkFree(double x, double y);
        Obstacle(bool* grid[PIXEL_HEIGHT]){
            for(int i=0;i<PIXEL_HEIGHT;i++){
                occupancy[i] = grid[i];
            }
        }
};

bool Obstacle::checkFree(int u, int v){
    return occupancy[v][u];
}
bool Obstacle::checkFree(double x, double y){
    // return true; // Ignore Obstacles;
    if ((x < 0) || (y < 0) || (x >= PIXEL_WIDTH/DENSITY) || (y>= PIXEL_HEIGHT/DENSITY)){
        return false;
    }
    return !occupancy[(int)(y*DENSITY)][(int)(x*DENSITY)];
}

class LineChecker {
    private:
        Obstacle* obst;
        bool checkLineHigh(int x0, int y0, int x1, int y1);
        bool checkLineLow(int x0, int y0, int x1, int y1);
    public:
        LineChecker(Obstacle* _obst){
            obst = _obst;
        }
        bool ObstacleFree(int x0, int y0, int x1, int y1);
};

bool LineChecker::checkLineHigh(int x0, int y0, int x1, int y1){
    bool isfree = true;
    int dx = x1 - x0;
    int dy = y1 - y0;
    int xi = 1;
    if (dx < 0) {
        xi = -1;
        dx = -dx;
    }
    int D = 2*dx - dy;
    int x = x0;
    int y;
    for(y=y0;y<=y1;y++){
        if (obst->checkFree(x, y)) {
            isfree = false;
            break;
        }
        if (D > 0){
            x = x + xi;
            D = D - 2*dy;
        }
        D = D + 2*dx;
    }
    return isfree;
}

bool LineChecker::checkLineLow(int x0, int y0, int x1, int y1){
    bool isfree = true;
    int dx = x1 - x0;
    int dy = y1 - y0;
    int yi = 1;
    if (dy < 0) {
        yi = -1;
        dy = -dy;
    }
    int D = 2*dy - dx;
    int y = y0;
    int x;
    for(x=x0;x<=x1;x++){
        if (obst->checkFree(x, y)) {
            isfree = false;
            break;
        }
        if (D > 0){
            y = y + yi;
            D = D - 2*dx;
        }
        D = D + 2*dy;
    }
    return isfree;
}

bool LineChecker::ObstacleFree(int x0, int y0, int x1, int y1){
    bool isfree = true;
    if (abs(y1 - y0) < abs(x1 - x0)){
        if (x0 > x1){
            isfree = checkLineLow(x1, y1, x0, y0);
        }
        else {
            isfree = checkLineLow(x0, y0, x1, y1);
        }
    }
    else {
        if (y0 > y1){
            isfree = checkLineHigh(x1, y1, x0, y0);
        }
        else {
            isfree = checkLineHigh(x0, y0, x1, y1);
        }
        
    }
    return isfree;
}

class Node {
    public:
        int i;
        double x, y, theta, cost;
        Node* parent;
        Node(){

        }
        Node(int _i, double _x, double _y, double _theta=0) {
            i = _i;
            x = _x;
            y = _y;
            theta = _theta;
            parent = nullptr;
        }
        double Distance(Node* node2);
        Node* Steer(Node* goal, Obstacle &obst, bool limit_length);
        void print(ofstream& file){
            if (parent == nullptr){
                file << "(" << i << ", None,(" << x << ", " << y <<")),\n";
            }
            else
            {
                file << "(" << i << ", "<< parent->i <<",(" << x << ", " << y <<")),\n";
            }
        }
};

double Node::Distance(Node* node2){
    double d2 = (x-node2->x)*(x-node2->x) + (y-node2->y)*(y-node2->y);
    double angle_diff = (theta-node2->theta)*(theta-node2->theta);
    angle_diff = min(angle_diff, (theta-node2->theta+2*M_PI)*(theta-node2->theta+2*M_PI));
    angle_diff = min(angle_diff, (theta-node2->theta-2*M_PI)*(theta-node2->theta-2*M_PI));
    d2 += angle_diff * ANGLE_DISTANCE_FACTOR;
    return sqrtf64(d2);
}

Node* Node::Steer(Node* goal, Obstacle &obst, bool limit_length = true){
    double dist = Distance(goal);
    Node* new_node;
    if (limit_length && dist > DT) {
        new_node = new Node(goal->i, x+(goal->x-x)*(DT/dist), y+(goal->y-y)*(DT/dist));
    } else {
        new_node = goal;
    }
    if (LineChecker(&obst).ObstacleFree(
        (int)(x*DENSITY), (int)(y*DENSITY), 
        (int)(new_node->x*DENSITY), (int)(new_node->y*DENSITY)
    )){
        return new_node;
    }
    return nullptr;
}

Node* global_goal;

void constrain(double &x){
    while (x > M_PI){
        x -= 2*M_PI;
    }
    while (x < -M_PI){
        x += 2*M_PI;
    }
}

class Control {
    public:
        double v;
        double omega;
        Control(double _v, double _omega){
            v = _v;
            omega = _omega;
        }
};

class State {
    public:
        double x;
        double y;
        double theta;

        double rho;
        double alpha;
        double phi;
        State(){

        }
        State(Node* node, Node* goal) {
            x = node->x;
            y = node->y;
            theta = node->theta;
            Cartesian2Polar(goal);
        } 
        void Cartesian2Polar(Node* goal){
            rho = (x-goal->x)*(x-goal->x)+(y-goal->y)*(y-goal->y);
            rho = sqrtf64(rho);
            alpha = atan2f64(goal->y-y, goal->x-x)-theta;
            phi = goal->theta-theta;
            constrain(alpha);
            constrain(phi);
        }
        void Polar2Cartesian(Node* goal){
            theta = goal->theta - phi;
            constrain(theta);
            x = goal->x - rho*cosf64(alpha+theta);
            y = goal->y - rho*sinf64(alpha+theta);
        }
        double CostMetric(State* s) {
            double c2 = (x-s->x)*(x-s->x) + (y-s->y)*(y-s->y);
            // Angles are chosed to be within pi/2 difference
            double angle_diff = 1-cosf64(theta-s->theta);
            c2 += angle_diff * angle_diff;
            return sqrtf64(c2);
        }
        double Distance(Node* s) {
            double d2 = (x-s->x)*(x-s->x) + (y-s->y)*(y-s->y);
            double angle_diff = (theta-s->theta)*(theta-s->theta);
            angle_diff = min(angle_diff, (theta-s->theta+M_PI)*(theta-s->theta+2*M_PI));
            angle_diff = min(angle_diff, (theta-s->theta-2*M_PI)*(theta-s->theta-2*M_PI));
            d2 += angle_diff * ANGLE_DISTANCE_FACTOR;
            return sqrtf64(d2);
        }
        State* step(Control control, Node* goal){
            State* newState = new State();
            // Using Euler Integration
            newState->rho = rho - (cosf64(alpha) * control.v)*DT;
            newState->alpha = alpha + (sin(alpha)*control.v/rho-control.omega)*DT;
            newState->phi = phi - control.omega*DT;
            newState->Polar2Cartesian(goal);
            return newState;
        }
};

Control* getControl(State s){
    
    double v = tanhf64(3.8*s.rho);
    double omega = 6 * s.alpha - s.phi;
    if (v > MAX_V) {
        v = MAX_V;
    }
    if (omega > MAX_OMEGA) {
        omega = MAX_OMEGA;
    }
    if (omega < -MAX_OMEGA) {
        omega = -MAX_OMEGA;
    }
    Control* control = new Control(v, omega);
    return control;
}

Node* SampleFree(int i, double x0, double x1, double y0, double y1){
    Node* node = new Node(i, randDouble(x0, x1), randDouble(y0, y1), randDouble(-M_PI, M_PI));
    return node;
}


class NHL_Node:public Node{
    public:
        vector<State*> states; // parent to this
        vector<Control*> controls; // parent to this
        NHL_Node():Node(){

        }
        NHL_Node(int _i, double _x, double _y, double _theta):Node(_i, _x, _y, _theta){

        }
        NHL_Node* Steer(Node* goal, Obstacle &obst, bool force_connect = false){
            NHL_Node* new_node = new NHL_Node();
            new_node->i = goal->i;
            State* s = new State(this, goal);
            new_node->states.push_back(s);
            bool connectable = false;
            for(int i=0;i<MAX_STEER_ITER;i++){
                if (s->Distance(goal) < GAMMA || s->Distance(global_goal) < GAMMA) {
                    // cout << "Success\n";
                    connectable = true;
                    break;
                }
                Control* c = getControl(*s);
                s = s->step(*c, goal);
                if(!obst.checkFree(s->x, s->y)){
                    // cout << "Collision\n";
                    return nullptr;
                }
                new_node->states.push_back(s);
                new_node->controls.push_back(c);
            }
            // reverse(new_node->states.begin(), new_node->states.end());
            // reverse(new_node->controls.begin(), new_node->controls.end());
            if (force_connect) {
                if (connectable) {
                    controls.push_back(getControl(*s));
                    State* goal_state = new State();
                    goal_state->x = goal->x;
                    goal_state->y = goal->y;
                    goal_state->theta = goal->theta;
                    new_node->states.push_back(goal_state);
                    // TODO Check obstacle for line here.
                    new_node->x = goal->x;
                    new_node->y = goal->y;
                    new_node->theta = goal->theta;
                    return new_node;
                } else {
                    return nullptr;
                }
            }
            new_node->x = s->x;
            new_node->y = s->y;
            new_node->theta = s->theta;
            return new_node;
        }
        void print(ofstream& file){
            if (parent == nullptr){
                file << "(" << i << ", None,(" << x << ", " << y <<"), None, None),\n";
            }
            else
            {
                file << "(" << i << ", "<< parent->i <<",(" << x << ", " << y <<"), (";
                for(int i=0; i<states.size(); i++){
                    file << states[i]->x << ", ";
                }
                file << "), (";
                for(int i=0; i<states.size(); i++){
                    file << states[i]->y << ", ";
                }
                file << ")),\n";
            }
        }
        double CostMetric(){
            double self_cost = 0;
            for(int i=0; i<states.size()-1; i++){
                self_cost += states[i]->CostMetric(states[i+1]);
            }
            return self_cost;
        }
};

NHL_Node* NHLSampleFree(int i, double x0, double x1, double y0, double y1){
    if(randDouble(0, 1) < SAMPLE_GOAL_RATE){
        return (NHL_Node*) global_goal;
    }
    NHL_Node* node = new NHL_Node(i, randDouble(x0, x1), randDouble(y0, y1), randDouble(-M_PI, M_PI));
    return node;
}

class Graph {
    public:
        vector<Node*> head;
        int n; // number of nodes
        Graph(Node* root){
            head.push_back(root);
            n = 1;
        }
        Node* Nearest(Node &node){
            Node* nearest = head[0];
            double minDist = node.Distance(nearest);
            for(int i=1;i<n;i++){
                double newDist = node.Distance(head[i]);
                if (newDist < minDist) {
                    nearest = head[i];
                    minDist = newDist;
                }
            }
            return nearest;
        }
        void FindNear(vector<Node*> &near, Node* node, Obstacle &obst){
            for(int i=0; i<n; i++){
                if (head[i]->Distance(node) < NEAR_RADIUS && 
                        head[i]->Steer(node, obst, false) != nullptr){
                    near.push_back(head[i]);
                }
            }
        }
        void printAll(ofstream& file){
            file << "\nnodes_list = (\n";
            for(int i=0; i<n;i++){
                head[i]->print(file);
            }
            file << ")\n";
        }
        void printAllNHL(ofstream& file){
            file << "\nnodes_list = (\n";
            for(int i=0; i<n;i++){
                ((NHL_Node*)head[i])->print(file);
            }
            file << ")\n";
        }
        void printAns(ofstream& file, Node* node = nullptr){
            file << "points = (\n";
            if (node == nullptr) {
                node = head[n-1];
            }
            while (node != nullptr) {
                file << "(" << node->x << ", " << node->y << "),\n";
                node  = node->parent;
            }
            file << ")\n";
        }
        void printAnsNHL(ofstream& file, Node* node = nullptr){
            file << "points = (\n";
            if (node == nullptr) {
                node = head[n-1];
            }
            while (node != nullptr) {
                ((NHL_Node*)node)->print(file);
                node  = node->parent;
            }
            file << ")\n";
        }
        void addNode(Node* node){
            head.push_back(node);
            n += 1;
        }
};


void RRT(Node init, Node goal, Obstacle obst, ofstream& file){
    Graph graph(&init);
    bool arrived = false;
    for(int i=1; i<N; i++) {
        Node* x_rand = SampleFree(graph.n, 0, (double)PIXEL_WIDTH/(double)DENSITY, 0, (double)PIXEL_HEIGHT/(double)DENSITY);
        Node* x_nearest = graph.Nearest(*x_rand);
        Node* x_new = x_nearest->Steer(x_rand, obst);
        if (x_new != nullptr) {
            x_new->parent = x_nearest;
            graph.addNode(x_new);
            if (x_new->Distance(&goal) < TERMINATE_DIST){
                arrived = true;
                break;
            }
        }
    }
    if (arrived){
        cout << "Goal achieved with exploring " << graph.n << " nodes\n";
        graph.printAns(file);
        graph.printAll(file);
    } else {
        cout << "Goal not achieved\n";
        graph.printAll(file);
    }
}

void NHL_RRT(NHL_Node init, NHL_Node goal, Obstacle obst, ofstream& file){
    Graph graph(&init);
    bool arrived = false;
    for(int i=1; i<N; i++) {
        NHL_Node* x_rand = NHLSampleFree(graph.n, 0, (double)PIXEL_WIDTH/(double)DENSITY, 0, (double)PIXEL_HEIGHT/(double)DENSITY);
        NHL_Node* x_nearest = (NHL_Node*) graph.Nearest(*x_rand);
        NHL_Node* x_new = x_nearest->Steer(x_rand, obst);
        if (x_new != nullptr) {
            x_new->parent = x_nearest;
            graph.addNode(x_new);
            if (x_new->Distance(&goal) < TERMINATE_DIST){
                arrived = true;
                break;
            }
        }
    }
    if (arrived){
        cout << "Goal achieved with exploring " << graph.n << " nodes\n";
        graph.printAnsNHL(file);
        graph.printAllNHL(file);
    }
    else {
        cout << "Goal not achieved\n";
        graph.printAllNHL(file);
    }
}

void RRT_star(Node init, Node goal, Obstacle obst, ofstream& file){
    init.cost = 0;
    Graph graph(&init);
    bool arrived = false;
    for(int i=1; i<N; i++) {
        Node* x_rand = SampleFree(graph.n, 0, (double)PIXEL_WIDTH/(double)DENSITY, 0, (double)PIXEL_HEIGHT/(double)DENSITY);
        Node* x_nearest = graph.Nearest(*x_rand);
        Node* x_new = x_nearest->Steer(x_rand, obst);
        if (x_new != nullptr) {
            Node* x_min = x_nearest;
            double c_min = x_nearest->cost + x_nearest->Distance(x_new);
            vector<Node*> X_near;
            graph.FindNear(X_near, x_new, obst);
            for(int j=0; j<X_near.size(); j++){
                Node* x_near = X_near[j];
                double new_cost = x_near->cost + x_near->Distance(x_new);
                if (new_cost < c_min) {
                    x_min = x_near;
                    c_min = new_cost;
                }
            }
            x_new->parent = x_min;
            x_new->cost = c_min;
            for(int j=0; j<X_near.size(); j++){
                Node* x_near = X_near[j];
                double t = x_new->cost + x_new->Distance(x_near);
                if (t < x_near->cost) {
                    x_near->parent = x_new;
                    x_near->cost = t;
                }
            }

            graph.addNode(x_new);
            if (x_new->Distance(&goal) < TERMINATE_DIST){
                arrived = true;
                break;
            }
        }
    }

    if (arrived){
        cout << "Goal achieved with exploring " << graph.n << " nodes\n";
        graph.printAns(file);
        graph.printAll(file);
    } else {
        cout << "Goal not achieved\n";
        graph.printAll(file);
    }
}

void NHL_RRT_star(NHL_Node init, NHL_Node goal, Obstacle obst, ofstream& file){
    init.cost = 0;
    Graph graph(&init);
    bool arrived = false;
    NHL_Node* shortest_end = nullptr;
    for(int i=1; i<N; i++) {
        NHL_Node* x_rand = NHLSampleFree(graph.n, 0, (double)PIXEL_WIDTH/(double)DENSITY, 0, (double)PIXEL_HEIGHT/(double)DENSITY);
        NHL_Node* x_nearest = (NHL_Node*) graph.Nearest(*x_rand);
        NHL_Node* x_new = x_nearest->Steer(x_rand, obst);
        if (x_new != nullptr) {
            NHL_Node* x_min = x_nearest;
            NHL_Node* x_min_new = nullptr;
            double c_min = x_nearest->cost + x_new->CostMetric();
            vector<NHL_Node*> X_near;
            vector<NHL_Node*> X_near_new;
            for(int j=0; j<graph.n; j++){
                if (graph.head[j]->Distance(x_new) < NEAR_RADIUS){ 
                    NHL_Node* x_new_near = x_new->Steer(graph.head[j], obst, true);
                    if(x_new_near == nullptr) {
                        continue;
                    }
                    X_near.push_back((NHL_Node*) graph.head[j]);
                    X_near_new.push_back(x_new_near);
                }
            }

            for(int j=0; j<X_near.size(); j++){
                NHL_Node* x_near = X_near[j];
                double new_cost = x_near->cost + X_near_new[j]->CostMetric();
                if (new_cost < c_min) {
                    x_min = x_near;
                    x_min_new = X_near_new[j];
                    c_min = new_cost;
                }
            }
            x_new->parent = x_min;
            x_new->cost = c_min;
            if (x_min_new != nullptr){
                // cout << "Init - Rewired " << x_new->i << endl;
                x_new->states = x_min_new->states;
                reverse(x_new->states.begin(), x_new->states.end());
                x_new->controls = x_min_new->controls;
                reverse(x_new->controls.begin(), x_new->controls.end());
            }

            for(int j=0; j<X_near.size(); j++){
                NHL_Node* x_near = X_near[j];
                double t = x_new->cost + X_near_new[j]->CostMetric();
                if (t < x_near->cost) {
                    // cout << "Post - Rewired " << x_near->i << endl;
                    x_near->parent = x_new;
                    x_near->cost = t;
                    x_near->states = X_near_new[j]->states;
                    x_near->controls = X_near_new[j]->controls;
                }
            }

            graph.addNode(x_new);
            if (x_new->Distance(&goal) < TERMINATE_DIST){
                arrived = true;
                if(shortest_end == nullptr) {
                    shortest_end = x_new;
                } else if (shortest_end->cost > x_new->cost) {
                    shortest_end = x_new;
                }
                // break;
            }
        }
    }

    if (arrived){
        cout << "Goal achieved with exploring " << graph.n << " nodes\n";
        graph.printAnsNHL(file, shortest_end);
        graph.printAllNHL(file);
    }else {
        cout << "Goal not achieved\n";
        graph.printAllNHL(file);
    }
}

int main(){
    srand((unsigned) time(0));
    // srand(0);
    ofstream file("result.py");
    NHL_Node init(0, 0.9, 0.9, M_PI_4+M_PI);
    NHL_Node goal(0, 0.1, 0.1, M_PI_4);
    global_goal = &goal;
    bool* obstacles[PIXEL_HEIGHT];
    for (int i=0;i<PIXEL_HEIGHT;i++){
        obstacles[i] = new bool[PIXEL_WIDTH];
        for (int j=0; j<PIXEL_WIDTH; j++){
            if(i>2 && i<7 && j>2 && j<7){
                obstacles[i][j] = true;
            } else {
                obstacles[i][j] = false;
            }
        }
    }
    file << "points = ()\n";

    Obstacle obst(obstacles);
    // RRT_star(init, goal, obst, file);
    // RRT(init, goal, obst, file);
    // NHL_RRT(init, goal, obst, file);
    int time = clock();
    NHL_RRT_star(init, goal, obst, file);
    int time_end = clock();
    cout << (time_end - time) / (double) CLOCKS_PER_SEC << endl;
    /*
    NHL_Node* new_node = init.Steer(&goal, obst, true);
    if (new_node == nullptr) {
        cout << "Routing Failed\n";
        file.close();
        return 0;
    
    }
    State s(&init, &goal);
    cout << s.x << " " << s.y << " " << s.theta << endl;
    // s.Polar2Cartesian(&goal);
    // cout << s.x << " " << s.y << endl;
    Control* control = getControl(s);
    cout << control->v << " " << control->omega << endl;
    s = *(s.step(*control, &goal));
    cout << s.x << " " << s.y << " " << s.theta << endl;
    cout << s.rho << " " << s.alpha << " " << s.phi << endl;
    new_node->parent = &init;
    file << "\nnodes_list = (\n";
    new_node->print(file);
    file << ")\n";
    */

    file.close();
    return 0;
}