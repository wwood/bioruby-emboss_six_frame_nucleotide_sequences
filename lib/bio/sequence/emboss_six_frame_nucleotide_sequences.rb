require 'bio'

module Bio
  class Sequence
    class NA
      # Translate this nucleotide sequence into a particular frame, but return
      # the nucleotide sequence of that translation, rather the amino acid
      # sequence. The frame returned aligns with the frames that are generated by
      # bioruby itself and the EMBOSS package program transeq.
      #
      #   Bio::Sequence::NA.new('ATGATG').transeq_nucleotide_sequence(1) => Bio::Sequence::NA.new('ATGATG')
      #   Bio::Sequence::NA.new('ATGATG').transeq_nucleotide_sequence(2) => Bio::Sequence::NA.new('TGATG')
      #   Bio::Sequence::NA.new('ATGATG').transeq_nucleotide_sequence(4) => Bio::Sequence::NA.new('CATCAT')
      def nucleotide_sequence_of_translation(frame)
        unless [-1,-2,-3,1,2,3,4,5,6].include?(frame) #error checking
          raise Exception, "unexpected frame for translation: `#{frame.inspect}'"
        end

        # Offset table for reverse frames. indexed by frame-4, then length%3
        offset_table = [[0,-2,-1],[-1,0,-2],[-2,-1,0]]

        # deal with the easy case of translating in the forward direction.
        if frame < 4
          return Bio::Sequence::NA.new(self[frame-1..length-1])
        end

        # translate negatives into positives for reverse sequences
        frame = 6 if frame == -3
        frame = 5 if frame == -2
        frame = 4 if frame == -1

        remainder = length%3
        offset = offset_table[remainder][frame-4]
        return Bio::Sequence::NA.new(self[0..length-1+offset].reverse_complement)
      end
    end
  end
end