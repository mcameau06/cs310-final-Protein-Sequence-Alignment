class NeedlemanWunsch:
    def __init__(self, seq1:str, seq2:str, match_score=1, mismatch_penalty=-1, gap_penalty=-2):
        self.seq1 = seq1
        self.seq2 = seq2
        self.match_score = match_score
        self.mismatch_penalty = mismatch_penalty
        self.gap_penalty = gap_penalty
        self.score_matrix = []

    def initialize(self):
        self.score_matrix = [[0 for i in range(len(self.seq2) + 1)] for j in range(len(self.seq1) + 1)]
        
        # add gap penalty to the first row and first column
        for i in range(len(self.seq1)+1):
          self.score_matrix[i][0] = i * self.gap_penalty
        
        for i in range(len(self.seq2)+1):
          self.score_matrix[0][i] = i * self.gap_penalty 

    def fill(self):
      diagonal_score = 0
      for i in range(1, len(self.seq1)+1):
        for j in range(1, len(self.seq2)+1):
          # check if there is a match
          if self.seq1[i-1] == self.seq2[j-1]:
            diagonal_score = self.match_score
          else:
            diagonal_score = self.mismatch_penalty

          # at each cell we compute what the optimal move will be
          self.score_matrix[i][j] = max(
              self.score_matrix[i-1][j-1] + diagonal_score, 
              self.score_matrix[i-1][j] + self.gap_penalty, 
              self.score_matrix[i][j-1] + self.gap_penalty)

    def propagate_backwards(self):
        aligned_1 = []
        aligned_2 = []
        diagonal_score = 0
        # start at the bottom right corner of the matrix
        i = len(self.seq1)
        j = len(self.seq2)
        while i > 0 or j > 0:

          if i == 0:
            aligned_1.append("-")
            aligned_2.append(self.seq2[j-1])
            j -= 1
          elif j == 0:
            aligned_1.append(self.seq1[i-1])
            aligned_2.append("-")
            i -= 1
          else:
            if self.seq1[i-1] == self.seq2[j-1]:
              diagonal_score = self.match_score
            else:
              diagonal_score = self.mismatch_penalty
            if self.score_matrix[i][j] == self.score_matrix[i-1][j-1] + diagonal_score:
              aligned_1.append(self.seq1[i-1])
              aligned_2.append(self.seq2[j-1])
              i , j = i-1, j-1
            elif self.score_matrix[i][j] == self.score_matrix[i-1][j] + self.gap_penalty:
              aligned_1.append(self.seq1[i-1])
              aligned_2.append('-')
              i , j = i-1, j
            elif self.score_matrix[i][j] == self.score_matrix[i][j-1] + self.gap_penalty:
              aligned_2.append(self.seq2[j-1])
              aligned_1.append('-')
              i , j = i, j-1

        return "".join(aligned_1[::-1]), "".join(aligned_2[::-1])

    def align(self):
        self.initialize()
        self.fill()
        return  self.propagate_backwards()

class SmithWaterman(NeedlemanWunsch):
  def __init__(self, seq1: str, seq2: str, match_score=1, mismatch_penalty=-1, gap_penalty=-2):
    super().__init__(seq1, seq2, match_score, mismatch_penalty, gap_penalty)

  def initialize(self):
        self.score_matrix = [[0 for i in range(len(self.seq2) + 1)] for j in range(len(self.seq1) + 1)]
        

  def fill(self):
      diagonal_score = 0
      for i in range(1, len(self.seq1)+1):
        for j in range(1, len(self.seq2)+1):
          # check if there is a match
          if self.seq1[i-1] == self.seq2[j-1]:
            diagonal_score = self.match_score
          else:
            diagonal_score = self.mismatch_penalty

          # at each cell we compute what the optimal move will be
          optimal_move = max(
              self.score_matrix[i-1][j-1] + diagonal_score, 
              self.score_matrix[i-1][j] + self.gap_penalty, 
              self.score_matrix[i][j-1] + self.gap_penalty, 0) # zero added to ensure opmimal move doesn't have negative cost
           
          self.score_matrix[i][j] = optimal_move

  def propagate_backwards(self):
      
      max_value = 0
      pos_i, pos_j = 0,0
      # find the maximum score for starting position
      for i in range(len(self.seq1)+1):
        for j in range(len(self.seq2) + 1):
          if self.score_matrix[i][j] > max_value:
            max_value = self.score_matrix[i][j]
            pos_i,pos_j = i,j
      aligned_1 = []
      aligned_2 = []
      diagonal_score = 0

      # start at the cell with the maximum score
      i,j = pos_i,pos_j

      while self.score_matrix[i][j] != 0:

          if i == 0:
            aligned_1.append("-")
            aligned_2.append(self.seq2[j-1])
            j -= 1
          elif j == 0:
            aligned_1.append(self.seq1[i-1])
            aligned_2.append("-")
            i -= 1

          else:
            if self.seq1[i-1] == self.seq2[j-1]:
              diagonal_score = self.match_score

            else:
              diagonal_score = self.mismatch_penalty
            if self.score_matrix[i][j] == self.score_matrix[i-1][j-1] + diagonal_score:
              aligned_1.append(self.seq1[i-1])
              aligned_2.append(self.seq2[j-1])

              i , j = i-1, j-1

            elif self.score_matrix[i][j] == self.score_matrix[i-1][j] + self.gap_penalty:
              aligned_1.append(self.seq1[i-1])
              aligned_2.append('-')

              i , j = i-1, j

            elif self.score_matrix[i][j] == self.score_matrix[i][j-1] + self.gap_penalty:
              aligned_2.append(self.seq2[j-1])
              aligned_1.append('-')

              i , j = i, j-1

      return "".join(aligned_1[::-1]), "".join(aligned_2[::-1])

  
