(ns bio-algs.proteomics
  (:require [clojure.set :refer [map-invert]]
            [clojure.string :refer [split]]
            [bio-algs.core :refer [read-file]]))

(def amino-acid-weights
  (zipmap
    ["G","A","S","P","V","T","C","I","L","N","D","Q","K","E","M","H","F","R","Y","W"]
    [ 57,71,87,97,99,101,103,113,113,114,115,128,128,129,131,137,147,156,163,186]))

(def weights-to-aa
  (map-invert amino-acid-weights))

(defn spec-graph
  "Find the graph for a given spectrum"
  [raw-spectrum]
  (let [spectrum (sort (conj raw-spectrum 0))
        n (count spectrum)]
    (reduce
      (fn [graph [w1 w2 aa]] (assoc-in graph [w1 w2] aa)) {}
      (for [i (range n) j (range (inc i) n)
            :let [w1 (nth spectrum i) w2 (nth spectrum j)
                  diff (- w2 w1) aa (weights-to-aa diff)]
            :when aa]
        [w1 w2 aa]))))

(defn draw-graph
  "Draws the graph back out"
  [graph]
  (let [sorted (sort (for [[start out] graph [end aa] out] [start end aa]))]
    (map (fn [[start end aa]] (str start "->" end ":" aa)) sorted)))

(defn all-paths
  "Finds all paths in the (acyclic) graph from one node to the other"
  [graph source sink]
  (defn- r-func
    [node]
    (if (= node sink) [""]
      (flatten
        (for [[next-node edge] (graph node)
              :let [tails (r-func next-node)]]
          (map #(str edge %) tails)))))
  (r-func source))

(defn weigh
  "weighs a given peptide"
  [peptide]
  (reduce + (map #(-> % str amino-acid-weights) peptide)))

(defn gen-spectrum
  [peptide]
  (let [n (count peptide)]
    (apply concat
           (for [i (range n)] (map weigh (split-at i peptide))))))

(defn explains-spectrum?
  [spectrum peptide]
  (= (sort (conj spectrum 0)) (sort (gen-spectrum peptide))))

(defn decode-spec-graph
  "Decodes the spectrum using the spectrum graph"
  [spectrum]
  (let [graph (spec-graph spectrum)
        paths (all-paths graph 0 (apply max spectrum))]
    (first (filter (partial explains-spectrum? spectrum) paths))))

(defn peptide-vector
  "Generates a peptide vector"
  ([peptide] (peptide-vector peptide amino-acid-weights))
  ([peptide aa-ws]
  (apply concat (for [aa peptide :let [w (aa-ws (str aa))]]
                  (reverse (cons 1 (repeat (dec w) 0)))))))

(defn from-pv
  "Generates the/a peptide that could generate the given peptide vector"
  ([pv] (from-pv pv weights-to-aa))
  ([pv ws-to-aa]
  (let [partitioned (partition 2 (partition-by identity pv))]
    (apply str (for [p partitioned] 
                 (ws-to-aa (inc (count (first p)))))))))

(defn dot-prod
  [v1 v2]
  (apply + (for [i (range (min (count v1) (count v2)))] (* (nth v1 i) (nth v2 i)))))

(defn score-peptide
  [peptide spectrum aa-ws]
  (let [n (count spectrum)
        pv (take n (peptide-vector peptide))
        match (apply str (take (apply + pv) peptide))]
    [(dot-prod spectrum pv) match]))

(defn best-peptide
  "Finds the best peptide matching a spectrum from a given protenome"
  ([spectrum proteome] (best-peptide spectrum proteome amino-acid-weights weights-to-aa))
  ([spectrum proteome aa-ws] (best-peptide spectrum proteome aa-ws (map-invert aa-ws)))
  ([spectrum proteome aa-ws ws-to-aa]
   (let [options (take-while not-empty (iterate rest proteome))
         scores (map #(score-peptide % spectrum aa-ws) options)
         best (apply max-key first scores)]
     (second best))))

(defn run
  [input]
  (let [in (read-file input)
        spectrum (map #(Integer/parseInt %) (split (first in) #"[ \t]"))
        graph (spec-graph spectrum)]
    (draw-graph graph)))
