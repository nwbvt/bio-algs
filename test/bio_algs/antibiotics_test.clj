(ns bio-algs.antibiotics-test
  (:require [clojure.test :refer :all]
            [bio-algs.antibiotics :refer :all]))

(deftest protein-translation
  (testing "That dna translates into proteins correctly"
    (is (= "MAMAPRTEINSTRING"
           (translate "AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA")))))

(deftest peptide-encoding
  (testing "Finding encodings for a peptide"
    (is (= ["ATGGCC" "GGCCAT" "ATGGCC"]
           (find-encoding "ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA" "MA")))))

(deftest theoretical-spectrum-test
  (testing "Finding the theoretical spectrum for a peptide"
    (is (= [0 113 114 128 129 227 242 242 257 355 356 370 371 484]
           (theoretical-spectrum "LEQN")))
    (is (= [0 113 114 128 129 242 242 257 370 371 484]
           (theoretical-spectrum "NQEL" true)))))

(deftest branch-and-bounds-test
  (testing "Using branch and bounds to perform sequencing"
    (is (= ["186-128-113" "186-113-128" "128-186-113" "128-113-186" "113-186-128" "113-128-186" ]
          (format-bb-seq [0 113 128 186 241 299 314 427]) ))))


(deftest cyclopeptide-scoring-test
  (testing "Scoring a peptide's experimental spectrum"
    (is (= 11 (score "NQEL" [0 99 113 114 128 227 257 299 355 356 370 371 484])))
    (is (== 8 (score "NQEL" [0 99 113 114 128 227 257 299 355 356 370 371 484] true)))
    ))

(deftest leaderboard-sequencing
  (testing "Sequencing via the leaderboard algorithm"
    (let [weights [0 71 113 129 147 200 218 260 313 331 347 389 460]]
     (is (= (score [113 147 71 129] weights)
           (score (lb-sequence 10 weights) weights))))))

(deftest convolution-test
  (testing "Testing the convolution of a spectrum"
    (is (= (sort [137 137 186 186 323 49])
           (sort (convolution [0 137 186 323]))))))

(deftest convolution-sequencing-test
  (testing "Testing using convolutions to sequence"
    (let [weights [57 57 71 99 129 137 170 186 194 208 228 265 285 299 307 323 356 364 394 422 493]]
      (is (= (score [99 71 137 57 72 57] weights))
          (score (convolution-seq 20 60 weights) weights)))))
