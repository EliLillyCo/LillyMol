class Accumulator
  attr_reader :n, :minval, :maxval, :sum
  def initialize
    @n = 0
    @minval = false
    @maxval = false
    @sum = false
    @sum2 = false
  end

  def reset
    @n = 0
    @minval = false
    @maxval = false
    @sum = false
    @sum2 = false
  end

  def extra(v)
    if (0 == self.n)
      @n = 1
      @minval = v
      @maxval = v
      @sum = v
      @sum2 = v * v
    else
      @n += 1
      @sum += v
      @sum2 += v * v
      if (v < self.minval)
        @minval = v
      elsif (v > @maxval)
        @maxval = v
      end
    end
  end

  def extra_n(v, n)
    if 0 == self.n
      @n = n
      @minval = v
      @maxval = v
      @sum = v
      @sum2 = v * v
    else
      @n += n
      @sum += n * v
      @sum2 += n * v * v
      if (v < self.minval)
        @minval = v
      elsif v > @maxval
        @maxval = v
      end
    end
  end

  def report(s)
    if (0 == self.n)
      s.print "NO data\n"
      return false
    end

    ave = @sum.to_f / @n.to_f

    s.print "#{@n} values between #{@minval} and #{@maxval} ave #{ave}\n"
  end

  def report(s, fmt)
    if (0 == self.n)
      s.print "NO data\n"
      return false
    end

    ave = @sum.to_f / @n.to_f

    printf_string = "#{@n} values between #{fmt} and #{fmt} ave #{fmt}\n"

    s.printf(printf_string, @minval, @maxval, ave)
  end

  def average
    return false if (0 == self.n)

    return @sum if (1 == self.n)

    return @sum.to_f / self.n.to_f
  end

  def variance
    return false if (@n < 2)

    ave = self.average

    return (@sum2 - @n * ave * ave) / (@n - 1)
  end

  def standard_deviation
    v = self.variance
    return false unless v
    return Math.sqrt(v)
  end
end
