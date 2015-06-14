(ns bio-algs.hmm
  (:require [bio-algs.core :refer :all]
            [clojure.string :refer [split trim]]))

(defn p-path
  [full-path trans-matrix]
  (loop [p 1/2 path full-path]
    (if (= 1 (count path)) p (recur (* p ((trans-matrix (first path)) (second path))) (rest path)))))

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
  (loop [current-state nil
         emissions (seq outcome)
         probability 1]
    (if (empty? emissions) probability
      (let [new-probabilities (zipmap states (map (fn [new-state]
                                                    (if (nil? current-state)
                                                      (* (/ 1 (count states))
                                                                         ((emission-matrix new-state) (first emissions)))
                                                      (apply + (map (fn [pre]
                                                                      (* (current-state pre)
                                                                         ((trans-matrix pre) new-state)
                                                                         ((emission-matrix new-state) (first emissions))))
                                                                    states))))
                                                  states))]
        (recur new-probabilities (rest emissions) (apply + (vals new-probabilities)))))))

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

(defn- missing?
  [sym]
  (= \- sym))

(defn- which-state?
  [advance? position [from-state sym]]
  (if advance?
    (if (missing? sym) (keyword (str \D (inc position)))
      (keyword (str \M (inc position))))
    (if (missing? sym) from-state
      (keyword (str \I position)))))

(defn- insert-state?
  [state]
  (.startsWith (str state) ":I"))

(defn- add-states
  [transitions old-states new-states advance?]
  (reduce (fn [trans [old-state new-state]]
            (if (or (insert-state? new-state) advance?)
              (assoc trans old-state (conj (or (trans old-state) [])
                                           new-state))
              trans))
          transitions
          (partition 2 (interleave old-states new-states))))

(defn- add-emissions
  [emissions states row advance?]
  (reduce (fn [ems [state em]]
            (if (and (not (missing? em))
                     (or advance? (insert-state? state)))
              (assoc ems state (conj (or (ems state) []) em))
              ems))
          emissions
          (partition 2 (interleave states row))))

(defn- compute-occurences
  [theta alignment]
  (let [n (count alignment)
        alignment-length (count (first alignment))]
    (loop [index 0 position 0
           transitions {}
           emissions {}
           last-states (repeat n :S)]
      (if (= alignment-length index) [(add-states transitions last-states (repeat n :E) true) emissions position]
        (let [col (map #(nth % index) alignment)
              advance? (> theta (/ (count (filter missing? col)) n))
              new-states (map (partial which-state? advance? position) (partition 2 (interleave last-states col)))
              new-transitions (add-states transitions last-states new-states advance?)
              new-emissions (add-emissions emissions new-states col advance?)]
          (recur (inc index) (if advance? (inc position) position) new-transitions new-emissions new-states)))))) 

(defn- probabilities
  [outcomes]
  (let [n (count outcomes)
        freqs (frequencies outcomes)]
    (zipmap (keys freqs) (map #(double (/ % n)) (vals freqs)))))

(defn- to-matrix
  [occurence-list]
  (loop [occurences occurence-list
         matrix {}]
    (if (empty? occurences) matrix
      (let [[from to] (first occurences)]
        (recur (rest occurences) (assoc matrix from (probabilities to)))))))

(defn profile-hmm
  [theta alignment]
  (let [[transitions emissions max-state] (compute-occurences theta alignment)
        transition-matrix (to-matrix transitions)
        emission-matrix (to-matrix emissions)]
    {:transition transition-matrix :emission emission-matrix :max-state max-state}))

(defn get-states
  [max-state]
  (flatten
    [[:S :I0]
     (for [i (range 1 (inc max-state))]
       (map #(keyword (str % i)) ["M" "D" "I"]))
    [:E]]))

(defn draw-row
  [matrix cols row]
  (apply str (name row) (map (fn [col]
                               (let [n (or ((or (matrix row) {}) col) 0)
                                     format-string (if (float? n) "\t%.3f" "\t%d")
                                     ]
                                 (format format-string n))) cols)))

(defn draw-profile
  [profile alphabet]
  (let [states (get-states (:max-state profile))]
    (flatten
      [(apply str (map #(format "\t%s" (name %)) states))
       (map (partial draw-row (profile :transition) states) states)
       "--------"
       (apply str (map #(format "\t%s" %) alphabet))
       (map (partial draw-row (profile :emission) alphabet) states)])))

(defn run
  [input-file]
  (let [in (vec (split (slurp input-file) #"\n"))
        theta (Double/parseDouble (first in))
        alphabet (map first (split (nth in 2) #"[ \t]"))
        alignment (drop 4 in)
        profile (profile-hmm theta alignment)]
    (draw-profile profile alphabet)))
