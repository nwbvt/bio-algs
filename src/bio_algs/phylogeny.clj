(ns bio-algs.phylogeny
  (:require [clojure.string :refer [split]]))

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
