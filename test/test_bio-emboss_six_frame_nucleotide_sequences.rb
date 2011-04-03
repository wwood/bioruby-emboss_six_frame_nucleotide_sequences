require 'helper'

class TestBioEmbossSixFrameNucleotideSequences < Test::Unit::TestCase
  na = Bio::Sequence::NA
  @@data_dir = File.join(File.dirname(__FILE__), 'data')

  should "test length divisible by 3 forwards" do
  # >9nt_1
  # MMM
  # >9nt_2
  # **X
  # >9nt_3
  # DDX
    a = na.new('ATGATGATG')
    assert_equal na.new('ATGATGATG').downcase, a.transeq_nucleotide_sequence(1)
    assert_equal na.new('TGATGATG').downcase, a.transeq_nucleotide_sequence(2)
    assert_equal na.new('GATGATG').downcase, a.transeq_nucleotide_sequence(3)
  end

  # >9nt_4
  # HHH
  # >9nt_5
  # SSX
  # >9nt_6
  # IIX
  should "test length divisible by 3 backwards" do
    a = na.new('ATGATGATG')
    assert_equal na.new('CATCATCAT').downcase, a.transeq_nucleotide_sequence(4)
    assert_equal na.new('TCATCAT').downcase, a.transeq_nucleotide_sequence(5)
    assert_equal na.new('ATCATCAT').downcase, a.transeq_nucleotide_sequence(6)
  end

  should "test length divisible by 3 remainder 1 forwards" do
  # >10nt_1
  # MMMX
  # >10nt_2
  # ***
  # >10nt_3
  # DDX
    a = na.new('ATGATGATGA')
    assert_equal na.new('ATGATGATGA').downcase, a.transeq_nucleotide_sequence(1)
    assert_equal na.new('TGATGATGA').downcase, a.transeq_nucleotide_sequence(2)
    assert_equal na.new('GATGATGA').downcase, a.transeq_nucleotide_sequence(3)
  end

  # >10nt_4
  # HHH
  # >10nt_5
  # SSSX
  # >10nt_6
  # IIX
  should "test length divisible by 3 remainder 1 backwards" do
    a = na.new('ATGATGATGA')
    assert_equal na.new('CATCATCAT').downcase, a.transeq_nucleotide_sequence(4)
    assert_equal na.new('TCATCATCAT').downcase, a.transeq_nucleotide_sequence(5)
    assert_equal na.new('ATCATCAT').downcase, a.transeq_nucleotide_sequence(6)
  end

  should "test length divisible by 3 remainder 2 forwards" do
  # >11nt_1
  # MMMX
  # >11nt_2
  # ***X
  # >11nt_3
  # DDD
    a = na.new('ATGATGATGAT')
    assert_equal na.new('ATGATGATGAT').downcase, a.transeq_nucleotide_sequence(1)
    assert_equal na.new('TGATGATGAT').downcase, a.transeq_nucleotide_sequence(2)
    assert_equal na.new('GATGATGAT').downcase, a.transeq_nucleotide_sequence(3)
  end

  # >11nt_4
  # HHH
  # >11nt_5
  # SSSX
  # >11nt_6
  # IIIX
  should "test length divisible by 3 remainder 2 backwards" do
    a = na.new('ATGATGATGAT')
    assert_equal na.new('CATCATCAT').downcase, a.transeq_nucleotide_sequence(4)
    assert_equal na.new('TCATCATCAT').downcase, a.transeq_nucleotide_sequence(5)
    assert_equal na.new('ATCATCATCAT').downcase, a.transeq_nucleotide_sequence(6)
  end

  # in test/data, 3 nucleotide sequences have been translated by transeq into
  # 9 different protein sequences. They should match the bioruby translations
  should "should align with the transeq translation" do
    nucleotide_sequences = {}
    protein_sequences = {}

    # Read in the files
    Bio::FlatFile.foreach(File.join(@@data_dir,'test.fa')) do |seq|
      nucleotide_sequences[seq.entry_id] = seq.seq
    end
    Bio::FlatFile.foreach(File.join(@@data_dir,'test_transeq_6frame.fa')) do |seq|
      protein_sequences[seq.entry_id] = seq.seq
    end

    # Make sure enough sequences are being tested
    assert_equal 3*6, protein_sequences.length

    # iterate them all and make sure they match
    protein_sequences.each do |pname, pseq|
      if matches = pname.match(/(.*)_([1-6])/)
        pseq.gsub!(/X/,'') #remove hanging Xs cos bioruby and transeq do that
        # differently
        assert_equal pseq, na.new(nucleotide_sequences[matches[1]]).transeq_nucleotide_sequence(matches[2].to_i).translate
      else
        raise
      end
    end
  end

end
