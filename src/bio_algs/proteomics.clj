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

(defn run
  [input]
  (let [in (read-file input)
        spectrum (map #(Integer/parseInt %) (split (first in) #"[ \t]"))
        graph (spec-graph spectrum)]
    (draw-graph graph)))
