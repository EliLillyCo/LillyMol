require_relative "group_args"
require "test/unit"

class TestGroupArgs < Test::Unit::TestCase

  def test_empty
    args = Array.new
    assert_nil(GfpMakeSupport.group_args(args))
  end

  def test_fails_non_option
    args = ["A", "AA"]
    result = GfpMakeSupport.group_args(args)
    assert_nil(result)
  end

  def test_fails_no_close_empty
    args = ["-MK", "-oB"]
    result = GfpMakeSupport.group_args(args)
    assert_nil(result)
  end

  def test_fails_no_close_not_empty
    args = ["-MK", "-oB", "a", "b", "-QQ"]
    result = GfpMakeSupport.group_args(args)
    assert_nil(result)
  end

  def test_unpaired_closing
    args = ["-AA", "-cB"]
    result = GfpMakeSupport.group_args(args)
    assert_nil(result)
  end

  def test_single_simple
    args = ["-QQ"]
    result = GfpMakeSupport.group_args(args)
    assert_equal(result.size, 1)
    assert_equal(result[0].option, "QQ")
    assert_nil(result[0].value)
  end

  def test_single_token_composite
    args = ["-oQQ", "a", "-cQQ"]
    result = GfpMakeSupport.group_args(args)
    assert_equal(result.size, 1)
    assert_equal(result[0].option, "QQ")
    assert_equal(result[0].value, "a")
  end

  def test_2_token_composite
    args = ["-oQQ", "a", "b", "-cQQ"]
    result = GfpMakeSupport.group_args(args)
    assert_equal(result.size, 1)
    assert_equal(result[0].option, "QQ")
    assert_equal(result[0].value, "a b")
  end

  def test_multi_token_composite
    args = ["-oQQ", "a", "b", "c", "d", "e", "-cQQ"]
    result = GfpMakeSupport.group_args(args)
    assert_equal(result.size, 1)
    assert_equal(result[0].option, "QQ")
    assert_equal(result[0].value, "a b c d e")
  end

  def test_single_and_composite_one
    args = ["-QQ", "-oQQ", "a", "-cQQ"]
    result = GfpMakeSupport.group_args(args)
    assert_equal(result.size, 2)
    assert_equal(result[0].option, "QQ")
    assert_nil(result[0].value)
    assert_equal(result[1].option, "QQ")
    assert_equal(result[1].value, "a")
  end

  def test_single_and_composite_one_rev
    args = ["-oQQ", "a", "-cQQ", "-QQ"]
    result = GfpMakeSupport.group_args(args)
    assert_equal(result.size, 2)
    assert_equal(result[0].option, "QQ")
    assert_equal(result[0].value, "a")
    assert_equal(result[1].option, "QQ")
    assert_nil(result[1].value)
  end

  def test_multiple_single
    args = "ABCDE".split("").map { |s| s.gsub(/^/, "-") }
    result = GfpMakeSupport.group_args(args)
    assert_equal(result.size, args.size)
    args.each_with_index do |arg, index|
      assert_equal(result[index].option, arg.gsub(/^-/, ""))
    end
  end

  def test_multiple_composite
    args = ["-oA", "1", "-cA", "-oA", "2", "3", "-cA", "-oB", "b", "b", "b", "-cB"]
    result = GfpMakeSupport.group_args(args)
    assert_equal(result.size, 3)
    assert_equal(result[0].option, "A")
    assert_equal(result[0].value, "1")
    assert_equal(result[1].option, "A")
    assert_equal(result[1].value, "2 3")
    assert_equal(result[2].option, "B")
    assert_equal(result[2].value, "b b b")
  end

  def test_mixture
    args = ["-oA", "11", "-cA", "-XX", "-Y", "-oA", "2", "3", "-cA", "-oB", "b", "b", "b", "-cB", "-CC"]
    result = GfpMakeSupport.group_args(args)
    assert_equal(result.size, 6)
    assert_equal(result[0].option, "A")
    assert_equal(result[0].value, "11")
    assert_equal(result[1].option, "XX")
    assert_nil(result[1].value)
    assert_equal(result[2].option, "Y")
    assert_nil(result[2].value)
    assert_equal(result[3].option, "A")
    assert_equal(result[3].value, "2 3")
    assert_equal(result[4].option, "B")
    assert_equal(result[4].value, "b b b")
    assert_equal(result[5].option, "CC")
    assert_nil(result[5].value)
  end
end

class TestGroupArgs < Test::Unit::TestCase
  def test_all_unique
    args = ["-A", "-B", "-oC", "q", "b", "-cC", "-D"]
    result = GfpMakeSupport.group_args(args)
    assert_equal(result.size, 4)
    unique = GfpMakeSupport.all_options_unique(result)
    assert_equal(unique.size, result.size)
  end
  def test_contains_simple_duplicate
    args = ["-A", "-A"]
    result = GfpMakeSupport.group_args(args)
    assert_equal(result.size, args.size)
    unique = GfpMakeSupport.all_options_unique(result)
    assert_nil(unique)
  end
  def test_contains_embedded_duplicate
    args = ["-A", "-B", "-A"]
    result = GfpMakeSupport.group_args(args)
    assert_equal(result.size, args.size)
    unique = GfpMakeSupport.all_options_unique(result)
    assert_nil(unique)
  end
  def test_contains_composite_duplicate
    args = ["-A", "-oB", "+", "-cB", "-A", "-B"]
    result = GfpMakeSupport.group_args(args)
    assert_equal(result.size, 4)
    unique = GfpMakeSupport.all_options_unique(result)
    assert_nil(unique)
  end
end
