(ns bio-algs.phylogeny
  (:require [clojure.string :refer [split]]
            [clojure.set :refer [intersection]]))

(defn- parse-text-edge
  "convert the textual format of an edge to an array"
  [text-edge]
  (let [[nodes weight] (split text-edge #":")
        [in-node out-node] (map keyword (split nodes #"->"))]
    [in-node out-node weight]))

(defn parse-graph
  "Parses a graph from the textual format to a map of nodes to an array of connected nodes and the weights of each"
  [text-edges]
  (let [edges (map parse-text-edge text-edges)
        edge-map (group-by first edges)
        nodes (keys edge-map)]
    (zipmap nodes
            (for [node nodes :let [edge-list (edge-map node)]]
              (zipmap (map second edge-list) (map #(Integer/parseInt (nth % 2)) edge-list))))))

(defn- is-leaf?
  "Returns true iff the given node is a leaf for the given graph"
  [graph node]
  (= 1 (count (graph node))))

(defn- int-val
  "Returns the integer value of the given keyword"
  [kw]
  (Integer/parseInt (name kw)))

(defn- distance
 "find the distance between two nodes in a graph"
 [graph node1 node2]
 (loop [nodes {node1 0}]
   (if (nil? (nodes node2))
     (if (empty? nodes) Double/POSITIVE_INFINITY
       (recur (apply merge (for [[node d] nodes
                                 :let [children (apply dissoc (graph node) (keys nodes))
                                       children-nodes (keys children)]]
                             (zipmap children-nodes (map #(+ d (children %)) children-nodes))))))
     (nodes node2))))

(defn weight-matrix
  "Find the distance matrix of all leaves to the other leaves"
  [graph]
  (let [leaves (sort-by int-val (filter (partial is-leaf? graph) (keys graph)))]
    (for [leaf1 leaves]
      (for [leaf2 leaves] (distance graph leaf1 leaf2)))))

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

(defn- adaptive-phylo
  "helper function for adaptive phylogeny (non tail) recursive algorithm"
  [matrix original-n]
  (let [n (dec (count matrix))
        n-k (keyword (str n))]
    (if (= 1 n) {:0 {:1 (get-in matrix [0 1])}
                 :1 {:0 (get-in matrix [1 0])}} 
      (let [limb (limb-length matrix n)
            nth-row (vec (map #(- % limb) (nth matrix n)))
            bald-partial (assoc matrix n nth-row)
            bald (vec (map-indexed (fn [i row] (assoc row n (nth-row i))) bald-partial))
            [i,j] (find-insertion-spot bald)
            [i-k, j-k] (map #(-> % str keyword) [i, j])
            new-matrix (subvec (vec (map #(subvec % 0 n) bald)) 0 n)
            t (adaptive-phylo new-matrix original-n)
            next-node (keyword (str (inc (+ original-n (- (count (keys t)) n)))))
            dist-to-i (- ((matrix i) n) limb)
            dist-to-j (- ((matrix j) n) limb)]
        (if (get-in t [i-k j-k])
            ;have to make the intervening node between two existing nodes
            (let [intervening next-node
                  old-i (t i-k) old-j (t j-k)
                  new-i (assoc (dissoc old-i j-k) intervening dist-to-i)
                  new-j (assoc (dissoc old-j i-k) intervening dist-to-j)]
              (assoc t i-k new-i j-k new-j intervening {i-k dist-to-i j-k dist-to-j n-k limb} n-k {intervening limb}))
            (if (and (t i-k) (t j-k))
              (let [intervening (first (intersection (set (keys (t i-k))) (set (keys (t j-k)))))]
                (assert intervening "Could not find intervening edge")
                (-> t (assoc-in [intervening n-k] limb)
                      (assoc-in [n-k intervening] limb)))
              (-> t (assoc-in [next-node n-k] limb)
                    (assoc-in [n-k next-node] limb)
                    (assoc-in [next-node i-k] dist-to-i)
                    (assoc-in [i-k next-node] dist-to-i)
                    (assoc-in [next-node j-k] dist-to-j)
                    (assoc-in [j-k next-node] dist-to-j)))))))) 

(defn make-tree
  "Make the tree matching the given distance matrix"
  [matrix]
  (adaptive-phylo matrix (count matrix)))
