(ns bio-algs.hmm)

(defn p-path
  [full-path trans-matrix]
  (loop [p 1/2 path full-path]
    (if (= 1 (count path)) p
      (recur (* p ((trans-matrix (first path)) (second path))) (rest path)))))
