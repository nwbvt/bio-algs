(ns bio-algs.core
  (:require [clojure.java.io :refer :all])
  (:require [clojure.string :refer [split trim join]]))

(defn write-result
  "Writes the solution in the expected format"
  [values]
  (with-open [outfile (writer "result.txt")]
    (binding [*out* outfile]
      (print (join "\n" values)))))

(defn write-1-result
  "Writes a single result in the expected format"
  [value]
  (write-result [value]))

(defn read-strands
  "Reads a list of \n deliminated dna strands"
  [input]
  (map trim (split input #"\n")))

(defn read-file
  "Reads strings from a file"
  [filename]
  (read-strands (slurp filename)))

(defn convert-lines
  "Change the results of a file read to a list of integers"
  [lines f]
  (vec (map (fn [line] (vec (map f (split line #"[ \t]")))) lines)))

(defn to-ints
  "Change the results of a file read to a list of integers"
  [lines]
  (convert-lines lines #(Integer/parseInt %)))

(defn to-floats
  "Change the results of a file read to a list of floating point numbers"
  [lines]
  (convert-lines lines #(Double/parseDouble %)))

(defn- parse-text-edge
  "convert the textual format of an edge to an array"
  [text-edge]
  (let [[nodes weight] (split text-edge #":")
        [in-node out-node] (map #(try (Integer/parseInt %) (catch NumberFormatException e %)) (split nodes #"->"))]
    [in-node out-node (if (nil? weight) nil (Double/parseDouble weight))]))

(defn parse-graph
  "Parses a graph from the textual format to a map of nodes to an array of connected nodes and the weights of each"
  [text-edges & {edge-parser :edge-parser}]
  (let [parse-text-edge (fn [text-edge]
                          (let [[nodes weight] (split text-edge #":")
                                [in-node out-node] (map #(try (Integer/parseInt %) (catch NumberFormatException e %))
                                                        (split nodes #"->"))]
                            [in-node out-node (if (nil? weight) nil ((or edge-parser  #(Double/parseDouble %)) weight))]))
        edges (map parse-text-edge text-edges)
        edge-map (group-by first edges)
        nodes (keys edge-map)]
    (zipmap nodes
            (for [node nodes :let [edge-list (edge-map node)]]
              (zipmap (map second edge-list) (map #(nth % 2) edge-list))))))


