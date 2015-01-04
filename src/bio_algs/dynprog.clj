(ns bio-algs.dynprog)

(defn find-change
  "given a money value and a list of coins, return the minimum number of coins needed to make that amount"
  [value coins]
  (loop [min-nums [0] v 1]
    (let [m (inc (apply min (for [c coins :when (<= c v)]
                              (min-nums (- v c)))))]
      (if (= value v) m
        (recur (conj min-nums m) (inc v))))))
