(ns bio-algs.phylogeny
  (:require [clojure.string :refer [split]]
            [clojure.set :refer [intersection]]))

(defn- parse-text-edge
  "convert the textual format of an edge to an array"
  [text-edge]
  (let [[nodes weight] (split text-edge #":")
        [in-node out-node] (map #(Integer/parseInt %) (split nodes #"->"))]
    [in-node out-node (Double/parseDouble weight)]))

(defn parse-graph
  "Parses a graph from the textual format to a map of nodes to an array of connected nodes and the weights of each"
  [text-edges]
  (let [edges (map parse-text-edge text-edges)
        edge-map (group-by first edges)
        nodes (keys edge-map)]
    (zipmap nodes
            (for [node nodes :let [edge-list (edge-map node)]]
              (zipmap (map second edge-list) (map #(nth % 2) edge-list))))))

(defn- is-leaf?
  "Returns true iff the given node is a leaf for the given graph"
  [graph node]
  (= 1 (count (graph node))))

(defn- walk-path
 "find the path and distance between two nodes"
 [graph node1 node2]
 (loop [nodes {node1 [[node1] 0.0]}]
   (if (nil? (nodes node2))
     (if (empty? nodes) [[] Double/POSITIVE_INFINITY]
       (recur (apply merge (for [[node [path d]] nodes
                                 :let [children (apply dissoc (graph node) (keys nodes))
                                       children-nodes (keys children)]]
                             (zipmap children-nodes (map #(vector (conj path %) (+ d (children %))) children-nodes))))))
     (nodes node2))))

(defn- distance
  "find the distance between two nodes in a graph"
  [graph node1 node2]
  (second (walk-path graph node1 node2)))

(defn- path
  "Find the path from one node to another"
  [graph node1 node2]
  (first (walk-path graph node1 node2)))

(defn weight-matrix
  "Find the distance matrix of all leaves to the other leaves"
  [graph]
  (let [leaves (sort (filter (partial is-leaf? graph) (keys graph)))]
    (vec (for [leaf1 leaves]
           (vec (for [leaf2 leaves] (distance graph leaf1 leaf2)))))))

(defn limb-length
  "Find the limb length of a node to its parent for graph with the given distance matrix"
  [matrix i]
  (let [n (count matrix)]
    (apply min (for [j (range n) k (range n) :when (and (not= i j) (not= i k))]
                 (/ (- (+ (-> matrix (nth i) (nth j)) (-> matrix (nth i) (nth k))) (-> matrix (nth j) (nth k))) 2)))))

(defn- find-insertion-spot
  "finds the insertion spot for the limb for the last leaf
  It will find two leaves i,j such that Dij = Din + Djn"
  [matrix]
  (let [n (dec (count matrix))]
    (first
      (for [i (range n) j (range n) 
            :when (= (get-in matrix [i j]) (+ (get-in matrix [i n]) (get-in matrix [j n])))]
        [i j]))))

(defn find-spot-in-path
  "finds the spot in the path the given distance away"
  [graph distance path]
  (loop [p path d distance]
    (let [a (first p) b (second p)
          part ((graph a) b)]
      (assert (not (empty? p)))
      (if (>= d part)
        (recur (rest p) (- d part))
        [[a b] d]))))

(defn- adaptive-phylo
  "helper function for adaptive phylogeny (non tail) recursive algorithm"
  [matrix original-n]
  (let [n (dec (count matrix))]
    (if (= 1 n) {0 {1 (get-in matrix [0 1])}
                 1 {0 (get-in matrix [1 0])}} 
      (let [limb (limb-length matrix n)
            nth-row (vec (map #(- % limb) (nth matrix n)))
            bald-partial (assoc matrix n nth-row)
            bald (vec (map-indexed (fn [i row] (assoc row n (nth-row i))) bald-partial))
            [node1 node2] (find-insertion-spot bald)
            new-matrix (subvec (vec (map #(subvec % 0 n) bald)) 0 n)
            t (adaptive-phylo new-matrix original-n)
            path (path t node1, node2)
            [[i j] d-to-i] (find-spot-in-path t (- ((matrix node1) n) limb) path)]
        (if (zero? d-to-i)
          ;just need to append to existing inner node i
          (-> t (assoc-in [i n] limb)
                (assoc-in [n i] limb))
          ;need to create a new innter node d away from i between i and j
          (let [next-node (+ original-n (- (count (keys t)) n))
                old-i (t i) old-j (t j) old-d ((t i) j)
                d-to-j (- old-d d-to-i)
                new-i (assoc (dissoc old-i j) next-node d-to-i)
                new-j (assoc (dissoc old-j i) next-node d-to-j)]
            (assoc t i new-i j new-j next-node {i d-to-i j d-to-j n limb} n {next-node limb}))))))) 

(defn make-tree
  "Make the tree matching the given distance matrix"
  [matrix]
  (adaptive-phylo matrix (count matrix)))

(defn format-graph
  "Formats the graph"
  [graph]
  (mapcat (fn [node] (map (fn [[k v]] (format "%d->%d:%f" node k v)) (graph node))) (sort (keys graph))))

(def cluster-dist
  "Find the distance between two clusters given a distance matrix"
  (memoize
    (fn [d-matrix c1 c2]
      (/ (apply + (for [m1 (:members c1) m2 (:members c2)]
                    ((d-matrix m1) m2)))
         (* (count (:members c1)) (count (:members c2)))))))

(defn- nearest-clusters
  "Find the clusters that are nearest given a distance matrix"
  [clusters d-matrix]
  (apply min-key first
         (for [c1 clusters c2 clusters :when (not= c1 c2)] [(cluster-dist d-matrix c1 c2) c1 c2])))

(defn upgma
  "Creates a ultrametric graph using unweighted pair group method with arithmetic mean"
  [matrix]
  (let [n (count matrix)]
   (loop [clusters (set (map #(hash-map :index % :members (list %)) (range n)))
          next-node n
          ages (vec (repeat n 0))
          t (zipmap (range n) (repeat n {}))]
     (if (= (count clusters) 1) t
       (let [[dist c1 c2] (nearest-clusters clusters matrix)
             age (/ dist 2)
             weight1 (double (- age (ages (:index c1))))
             weight2 (double (- age (ages (:index c2))))]
         (recur (conj (disj clusters c1 c2)
                      {:index next-node :members (concat (:members c1) (:members c2))})
                (inc next-node)
                (conj ages age)
                (-> t
                    (assoc-in [(:index c1) next-node] weight1)
                    (assoc-in [next-node (:index c1)] weight1)
                    (assoc-in [(:index c2) next-node] weight2)
                    (assoc-in [next-node (:index c2)] weight2))))))))
