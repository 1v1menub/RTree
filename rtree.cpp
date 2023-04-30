//#include <SDL2/SDL.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include <limits>
#include <cmath>
#include <cstdlib>
#include <ctime>

const int M = 3;
const int m = 2; // m <= (M/2)techo
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

    // Calcula el area necesaria para que el punto se encuentre contenido en el MBB
    double expansion_needed(const Point& p) const {
        double expansion = 0;
        for (int d = 0; d < NUM_DIMENSIONS; d++) {
            expansion += std::max(0.0, p.coords[d] - upper.coords[d]);
            expansion += std::max(0.0, lower.coords[d] - p.coords[d]);
        }
        return expansion;
    }

    // Permite expandir el MBB para que contenga al punto
    void expand(const Point& p) {
        for (int d = 0; d < NUM_DIMENSIONS; d++) {
            lower.coords[d] = std::min(lower.coords[d], p.coords[d]);
            upper.coords[d] = std::max(upper.coords[d], p.coords[d]);
        }
    }

    // Calcula el area del MBB
    double area() const {
        double area = 1;
        for (int d = 0; d < NUM_DIMENSIONS; d++) {
            area *= upper.coords[d] - lower.coords[d];
        }
        return area;
    }
};

// Distancia al cuadrado en todas las dimensiones => sirve para luego sacar pitagoras
static double point_distance_squared(const Point& p1, const Point& p2) {
    double distance_squared = 0;
    for (int d = 0; d < NUM_DIMENSIONS; d++) {
        double diff = p1.coords[d] - p2.coords[d];
        distance_squared += diff * diff;
    }
    return distance_squared;
}

// Distancia entre dos puntos => pitagoras
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
    std::shared_ptr<RTreeNode> parent_node() const { return parent; }
    void insert(const Point& p);
    
    static SplitType split_type;
    
private:
    bool is_leaf_;
    MBB mbb_;
    std::vector<Point> points_;
    std::vector<std::shared_ptr<RTreeNode>> children_;
    std::shared_ptr<RTreeNode> parent;

    void split();
    void linear_split();
    void quadratic_split();
    void brownie_split();
    std::shared_ptr<RTreeNode> choose_subtree(const Point& p) const;
    double compute_overlap(const MBB& mbb1, const MBB& mbb2) const;
};

RTreeNode::RTreeNode(bool is_leaf) : is_leaf_(is_leaf) {
    if (is_leaf) {
        mbb_ = MBB();
    }
    parent = nullptr;
}

void RTreeNode::insert(const Point& p) {
    if (is_leaf_) {
        // If the current node is a leaf, insert the point
        points_.push_back(p);
        mbb_.expand(p);

        // Check if the node overflows
        if (points_.size() > M) {
            // Split the node
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

    // Recorremos todos los hijos para escoger el que tenga menor expansion
    for (const std::shared_ptr<RTreeNode>& child : children_) {
        double expansion = child->mbb().expansion_needed(p);
        if (expansion < min_expansion) {
            min_expansion = expansion;
            best_child = child;
        } else if (expansion == min_expansion) {
            // En caso de empate, escogemos el hijo con menor area
            double area = child->mbb().area();
            if (area < best_child->mbb().area()) {
                best_child = child;
            }
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
    // Referencia a este mismo nodo
    std::shared_ptr<RTreeNode> thisptr = std::make_shared<RTreeNode>(this);

    if (this->is_leaf_) {
        // Crear un nuevo nodo hermano
        std::shared_ptr<RTreeNode> brother = std::make_shared<RTreeNode>(true);
        
        // Verificamos si existe padre
        if(this->parent) {
            parent = this->parent;
        } else{
            // Si no existe padre, creamos uno
            this->parent = std::make_shared<RTreeNode>(false);
        }

        // Asignamos el padre al nodo hermano
        brother->parent = this->parent;

        // Escoger los dos puntos más alejados entre sí
        double max_distance = 0;
        Point p1, p2;
        for (size_t i = 0; i < points_.size(); i++) {
            for (size_t j = i + 1; j < points_.size(); j++) {
                double distance = point_distance(points_[i], points_[j]);
                if (distance > max_distance) {
                    max_distance = distance;
                    p1 = points_[i];
                    p2 = points_[j];
                }
            }
        }

        // Creamos un vector de puntos temporal
        std::vector<Point> tempPoints = this->points_;
        // Limpiamos todos los puntos del nodo actual
        this->points_.clear();

        // Asignamos los puntos mas alejados a cada nodo
        this->insert(p1);
        brother->insert(p2);

        // Asignamos los hijos al nuevo padre
        parent->children_.push_back(thisptr);
        parent->children_.push_back(brother);

        // Agregar los puntos restantes al nodo cuya MBB tenga menor expansión
        for (const Point& p : tempPoints) {
            // Si el punto es igual a alguno de los puntos mas alejados, lo saltamos
            if (p.coords[0] == p1.coords[0] && p.coords[1] == p1.coords[1]) {
                continue;
            }
            if (p.coords[0] == p2.coords[0] && p.coords[1] == p2.coords[1]) {
                continue;
            }
            
            double expansion1 = this->mbb().expansion_needed(p);
            double expansion2 = brother->mbb().expansion_needed(p);

            if (expansion1 < expansion2) {
                this->insert(p);
            } else if (expansion2 < expansion1) {
                brother->insert(p);
            } else {
                // En caso de empate, escogemos el hijo con menor area de expansion necesaria
                double area1 = this->mbb().area();
                double area2 = brother->mbb().area();
                if (area1 < area2) {
                    this->insert(p);
                } else {
                    brother->insert(p);
                }
            }
        }
    }
}

void RTreeNode::quadratic_split(){
    // Se crean dos nodos hijos
    std::shared_ptr<RTreeNode> node1 = std::make_shared<RTreeNode>(is_leaf_);
    std::shared_ptr<RTreeNode> node2 = std::make_shared<RTreeNode>(is_leaf_);
    
    // Elegir los puntos que generan la mayor area entre sí
    double max_area = 0;
    Point p1, p2;
    for (size_t i = 0; i < points_.size(); i++) {
        for (size_t j = i + 1; j < points_.size(); j++) {
            double area = point_distance(points_[i], points_[j]);
            if (area > max_area) {
                max_area = area;
                p1 = points_[i];
                p2 = points_[j];
            }
        }
    }

    node1->insert(p1);
    node2->insert(p2);

    // Agregar los puntos restantes al nodo cuya MBB tenga menor area

    children_.push_back(node1);
    children_.push_back(node2);

    // Actualizar el MBB del nodo actual
    mbb_ = MBB();

    for (const Point& p : points_) {
        mbb_.expand(p);
    }

    node1->mbb_ = MBB();
    for (const Point& p : node1->points()) {
        node1->mbb_.expand(p);
    }

    // Limpiar los puntos del nodo actual
    points_.clear();

    // Actualizar el nodo actual como nodo interno
    is_leaf_ = false;

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
    //void draw(const RTreeNode* node, SDL_Renderer* renderer) const;
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
}
*/
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
/*
    // Initialize SDL
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