(ns bio-algs.hmm
  (:require [bio-algs.core :refer :all]
            [clojure.string :refer [split trim]]))

(defn p-path
  [full-path trans-matrix]
  (loop [p 1/2 path full-path]
    (if (= 1 (count path)) p
      (recur (* p ((trans-matrix (first path)) (second path))) (rest path)))))

(defn p-outcome
  [outcome path emission-matrix]
  (assert (= (count outcome) (count path)))
  (reduce * (map (fn [[emission state]] ((emission-matrix state) emission))
                 (partition 2 (interleave outcome path)))))

(defn viterbi
  "Uses the viterbi algorithm to find the optimal path"
  [outcome states trans-matrix emission-matrix]
  (let [cost (memoize (fn [state emission prior-state]
                        (* (:cost prior-state)
                           ((trans-matrix (:state prior-state)) state)
                           ((emission-matrix state) emission))))
        next-states (fn [emission previous-states]
                      (zipmap states (map (fn [state]
                                            (let [prior (previous-states
                                                          (apply max-key
                                                                 #(cost state emission
                                                                        (previous-states %))
                                                                 states))]
                                              {:cost (cost state emission prior)
                                               :state state
                                               :path (conj (:path prior) state)})) states)))]
    (loop [current (zipmap states
                           (map (fn [state] {:cost 1 :state state :path []}) states))
           emissions (seq outcome)]
      (if (empty? emissions)
        (let [best (apply max-key #((current %) :cost) states)]
          (apply str ((current best) :path)))
        (recur (next-states (first emissions) current) (rest emissions))))))

(defn optimal-path
  "Uses the viterbi algorithm to find the optimal path"
  [outcome states trans-matrix emission-matrix]
  (viterbi outcome states trans-matrix emission-matrix))

(defn p-hmm-outcome
  "Uses the viterbi algorithm to find the probability of an outcome given a hmm"
  [outcome states trans-matrix emission-matrix]
  nil)

(defn make-matrix
  [matrix-strings]
  (let [states (map first (split (trim (first matrix-strings)) #"[ \t]"))]
    (loop [rows (rest matrix-strings)
           matrix {}]
      (if (empty? rows) matrix
        (let [row (split (trim (first rows)) #"[ \t]")
              row-sym (first (first row))
              values (map #(Double/parseDouble %) (rest row))]
          (recur (rest rows) (assoc matrix row-sym (zipmap states values))))))))

(defn run
  [input-file]
  (let [in (vec (split (slurp input-file) #"\n"))
        outcome (first in)
        states (map first (split (in 4) #"[ \t]"))
        matrix-lines (drop 8 in)
        end-trans (+ 7 (count states))
        trans-matrix (make-matrix (subvec in 6 end-trans))
        emission-matrix (make-matrix (drop (inc end-trans) in))
        ]
    (viterbi outcome states trans-matrix emission-matrix)))
