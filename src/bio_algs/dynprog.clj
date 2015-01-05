(ns bio-algs.dynprog
  (:require [clojure.string :refer [split]]))

(defn find-change
  "given a money value and a list of coins, return the minimum number of coins needed to make that amount"
  [value coins]
  (loop [min-nums [0] v 1]
    (let [m (inc (apply min (for [c coins :when (<= c v)]
                              (min-nums (- v c)))))]
      (if (= value v) m
        (recur (conj min-nums m) (inc v))))))

(defn manhattan
  "Finds the longest distance in the manhatten tourist problem"
  [n m down right]
  (let [s (atom {})
        helper (fn [s i j]
                 (let [d-val (if (zero? i) 0 (+ (nth (nth down (dec i)) j) (s [(dec i) j])))
                      r-val (if (zero? j) 0 (+ (nth (nth right i) (dec j)) (s [i, (dec j)])))
                      v (max d-val r-val)]
                  (assoc s [i,j] v)))]
    (doseq [i (range 0 (inc n)) j (range 0 (inc m))]
      (swap! s helper i j))
    (@s [n, m])))

(defn read-in-manhattan
  [file]
  (let [make-ints (fn [strings] (map #(Integer/parseInt %) strings))
        lines (map #(make-ints (split % #" ")) (split (slurp file) #"\n"))
        [n, m] (first lines)
        down (take n (drop 1 lines))
        right (drop (inc n) lines)]
    {:n n :m m :down down :right right}))
