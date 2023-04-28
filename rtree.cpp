/* #include <SDL2/SDL.h> */
#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <ctime>

const int M = 3;
const int m = 2; // m <= M/2
const double EPS = 1e-9;
const int NUM_DIMENSIONS = 2;

enum SplitType {
    LINEAR_SPLIT,
    QUADRATIC_SPLIT,
    BROWNIE_SPLIT
};

struct Point {
    double coords[NUM_DIMENSIONS];
    Point() {}
    Point(double x, double y) {
        coords[0] = x;
        coords[1] = y;
    }
};

struct MBB {
    Point lower, upper;
    MBB(){
        lower = Point(std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
        upper = Point(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());
    }
    MBB(const Point& p1, const Point& p2) {
        for (int d = 0; d < NUM_DIMENSIONS; d++) {
            lower.coords[d] = std::min(p1.coords[d], p2.coords[d]);
            upper.coords[d] = std::max(p1.coords[d], p2.coords[d]);
        }
    }

    double expansion_needed(const Point& p) const {
        double expansion = 0;
        for (int d = 0; d < NUM_DIMENSIONS; d++) {
            expansion += std::max(0.0, p.coords[d] - upper.coords[d]);
            expansion += std::max(0.0, lower.coords[d] - p.coords[d]);
        }
        return expansion;
    }

    void expand(const Point& p) {
        for (int d = 0; d < NUM_DIMENSIONS; d++) {
            lower.coords[d] = std::min(lower.coords[d], p.coords[d]);
            upper.coords[d] = std::max(upper.coords[d], p.coords[d]);
        }
    }

    double area() const {
        double area = 1;
        for (int d = 0; d < NUM_DIMENSIONS; d++) {
            area *= upper.coords[d] - lower.coords[d];
        }
        return area;
    }
};


static double point_distance_squared(const Point& p1, const Point& p2) {
    double distance_squared = 0;
    for (int d = 0; d < NUM_DIMENSIONS; d++) {
        double diff = p1.coords[d] - p2.coords[d];
        distance_squared += diff * diff;
    }
    return distance_squared;
}

static double point_distance(const Point& p1, const Point& p2) {
    return std::sqrt(point_distance_squared(p1, p2));
}

class RTreeNode {
public:
    RTreeNode(bool is_leaf = true);
    
    bool is_leaf() const { return is_leaf_; }
    const MBB& mbb() const { return mbb_; }
    const std::vector<Point>& points() const { return points_; }
    const std::vector<std::shared_ptr<RTreeNode>>& children() const { return children_; }

    void insert(const Point& p);
    
    static SplitType split_type;
    
private:
    bool is_leaf_;
    MBB mbb_;
    std::vector<Point> points_;
    std::vector<std::shared_ptr<RTreeNode>> children_;
    RTreeNode* parent_;

    void split();
    void linear_split();
    void quadratic_split();
    void brownie_split();
    std::shared_ptr<RTreeNode> choose_subtree(const Point& p) const;
    /* double compute_overlap(const MBB& mbb1, const MBB& mbb2) const; */
};

RTreeNode::RTreeNode(bool is_leaf) : is_leaf_(is_leaf) {
    if (is_leaf) {
        mbb_ = MBB();
    }
}

void RTreeNode::insert(const Point& p) {
    if (is_leaf_) {
        points_.push_back(p);
        mbb_.expand(p);

        // Check if the node overflows
        if (points_.size() > M) {
            split();
        }
    } else {
        // If the current node is not a leaf, choose the best subtree to insert the point
        std::shared_ptr<RTreeNode> subtree = choose_subtree(p);
        subtree->insert(p);
        mbb_.expand(p);
    }
}

std::shared_ptr<RTreeNode> RTreeNode::choose_subtree(const Point& p) const {
    std::shared_ptr<RTreeNode> best_child;
    double min_expansion = std::numeric_limits<double>::max();
    
    for(auto i : this->children_) {
        double expand = i->mbb_.expansion_needed(p);
        if(expand < min_expansion) {
            min_expansion = expand;
            best_child = i;
        }
        else if(expand == min_expansion){

        }
    }
    return best_child;
}

void RTreeNode::split() {
    switch (split_type) {
        case LINEAR_SPLIT:
            linear_split();
            break;
        case QUADRATIC_SPLIT:
            quadratic_split();
            break;
        case BROWNIE_SPLIT:
            brownie_split();
            break;
        default:
            throw std::runtime_error("Unknown split type");
    }
}

SplitType RTreeNode::split_type;

void RTreeNode::linear_split() {
    // -------
    // Tu codigo aqui
    // -------

    std::shared_ptr<RTreeNode> thisptr = std::make_shared<RTreeNode>(this);
    if(this->is_leaf_) {
        Point initialp = this->points_[0];
        double bestdistance = -1000;
        int firstp;
        int lastp;
        Point p1;
        Point p2;
        for(int i = 0; i < this->points_.size(); i++) {
            double dist = point_distance(initialp, this->points_[i]);
            if(dist > bestdistance) {
                firstp = i;
                bestdistance = dist;
            }
        }
        p1 = points_[firstp];
        this->points_.erase(points_.begin() + firstp);
        bestdistance = 0;
        for(int i = 0; i < this->points_.size(); i++) {
            double dist = point_distance(this->points_[firstp], this->points_[i]);
            if(dist > bestdistance) {
                lastp = i;
                bestdistance = dist;
            }
        }
        p2 = points_[lastp];
        this->points_.erase(points_.begin() + lastp);


        std::vector<Point> puntos = this->points_;
        this->points_.clear();

        RTreeNode* parent;
        if(this->parent_) {
            parent = this->parent_;
        }
        std::shared_ptr<RTreeNode> sibling;

        this->insert(p1);
        sibling->insert(p2);

        this->parent_ = parent;
        sibling->parent_ = parent;

        parent->children_.push_back(thisptr);
        parent->children_.push_back(sibling);


        double ex1;
        double ex2;
        for(auto i : puntos) {
            ex1 = this->mbb_.expansion_needed(i);
            ex2 = sibling->mbb_.expansion_needed(i);
            if(ex1 < ex2) {
                this->insert(i);
            }
            else if(ex1 > ex2) {
                sibling->insert(i);
            }
            else {
                if(this->points_.size() < sibling->points_.size()) {
                    this->insert(i);
                }
                else {
                    sibling->insert(i);
                }
            }
            this->is_leaf_ = false;
        }
    }
    else {
        // IMEPLEMTNAR SPLIT PARA NODOS INTERNOS 
    }
}




void RTreeNode::quadratic_split(){
    // -------
    // Tu codigo aqui
    // -------
}
void RTreeNode::brownie_split(){
    // -------
    // Tu codigo aqui
    // -------
}


class RTree {
public:
    RTree(SplitType split_type = LINEAR_SPLIT);
    const RTreeNode* get_root() const { return root_.get(); }
    /* void draw(const RTreeNode* node, SDL_Renderer* renderer) const; */
    void print_ascii() const;
    void insert(const Point& p);

private:
    std::shared_ptr<RTreeNode> root_;
    void print_ascii_node(const RTreeNode* node, int depth = 0) const;
};

RTree::RTree(SplitType split_type) {
    RTreeNode::split_type = split_type;
    root_ = std::make_shared<RTreeNode>(true);
}

void RTree::insert(const Point& p) {
    root_->insert(p);
}
/* 
void RTree::draw(const RTreeNode* node, SDL_Renderer* renderer) const {
    if (node->is_leaf()) {
        SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
    } else {
        SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255);
    }

    SDL_Rect rect;
    rect.x = static_cast<int>(node->mbb().lower.coords[0]);
    rect.y = static_cast<int>(node->mbb().lower.coords[1]);
    rect.w = static_cast<int>(node->mbb().upper.coords[0] - node->mbb().lower.coords[0]);
    rect.h = static_cast<int>(node->mbb().upper.coords[1] - node->mbb().lower.coords[1]);

    SDL_RenderDrawRect(renderer, &rect);

    // Draw the points in the leaf nodes
    if (node->is_leaf()) {
        SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);

        for (const Point& p : node->points()) {
            // Draw a small rectangle (2x2 pixels) for each point
            SDL_Rect point_rect;
            point_rect.x = static_cast<int>(p.coords[0]) - 1;
            point_rect.y = static_cast<int>(p.coords[1]) - 1;
            point_rect.w = 2;
            point_rect.h = 2;
            SDL_RenderFillRect(renderer, &point_rect);
        }
    }
    for (const auto& child : node->children()) {
        draw(child.get(), renderer);
    }
} */

void RTree::print_ascii() const {
    print_ascii_node(root_.get());
}

void RTree::print_ascii_node(const RTreeNode* node, int depth) const {
    if (node == nullptr) {
        return;
    }

    std::string indentation(depth * 2, ' ');

    std::cout << indentation << "MBB: ("
              << node->mbb().lower.coords[0] << ", " << node->mbb().lower.coords[1] << "), ("
              << node->mbb().upper.coords[0] << ", " << node->mbb().upper.coords[1] << ")\n";

    if (node->is_leaf()) {
        std::cout << indentation << "Points:\n";
        for (const Point& p : node->points()) {
            std::cout << indentation << "  (" << p.coords[0] << ", " << p.coords[1] << ")\n";
        }
    } else {
        std::cout << indentation << "Children:\n";
        for (const auto& child : node->children()) {
            print_ascii_node(child.get(), depth + 1);
        }
    }
}



int main() {
    srand(time(NULL));
    RTree rtree(LINEAR_SPLIT);

    // Insert points
    rtree.insert(Point(  0,   0));
    rtree.insert(Point( 10,  10));
    rtree.insert(Point( 20,  20));
    rtree.insert(Point( 30,  30));
    rtree.insert(Point( 40,  40));
    rtree.insert(Point( 50,  50));
    rtree.insert(Point( 60,  60));
    rtree.insert(Point( 70,  70));
    rtree.insert(Point( 80,  80));
    rtree.insert(Point( 90,  90));
    rtree.insert(Point(100, 100));
    
    printf("Done\n");

    rtree.print_ascii();

   /*  // Initialize SDL
    if (SDL_Init(SDL_INIT_VIDEO) < 0) {
        std::cerr << "SDL initialization failed: " << SDL_GetError() << std::endl;
        return 1;
    }

    int window_width = 800;
    int window_height = 600;

    // Create a window
    SDL_Window* window = SDL_CreateWindow("R-Tree Visualization", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                          window_width, window_height, SDL_WINDOW_SHOWN);
    if (!window) {
        std::cerr << "Window creation failed: " << SDL_GetError() << std::endl;
        SDL_Quit();
        return 1;
    }

    // Create a renderer
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    if (!renderer) {
        std::cerr << "Renderer creation failed: " << SDL_GetError() << std::endl;
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    // R-Tree
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
    SDL_RenderClear(renderer);
    rtree.draw(rtree.get_root(), renderer);
    SDL_RenderPresent(renderer);

    // Wait 
    bool quit = false;
    SDL_Event event;
    while (!quit) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                quit = true;
            }
            if (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_ESCAPE) {
                quit = true;
            }
        }

        SDL_Delay(10); // Add a small delay to reduce CPU usage
    }

    // Clean up and exit
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
 */
    return 0;
}
