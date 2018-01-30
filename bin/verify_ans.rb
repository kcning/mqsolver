#!/usr/bin/env ruby
# verify a solution for an MQ problem
# usage: <challenge file> [<value for x_1> <value for x_2>, ... <value for x_n>]

sys_file, sol = ARGV[0..1]
sol = sol.delete('[ ]').split(',').map(&:to_i)
var_num = sol.size

sol.each do |s| 
  abort "wrong solution value: #{s}" if (s != 0 and s != 1)
end

mono_num = var_num * (var_num-1) / 2 + 2 * var_num + 1

# read the MQ system
File.open(sys_file, 'r') do |f|
  start = false
  counter = 0
  f.each_line do |line|
    if line =~ /\*+/
      start = true
      next
    end

    if start
      # parse the line, eq in grlex order
      eq = line.delete(';').split(' ').map(&:to_i)
      abort "wrong number of monomials in eq#{counter}: #{eq.size}" if eq.size != mono_num

      # evaluate the eq with the solution
      res = eq.pop # constant

      # linear terms
      var_num.times do |i|
        res += sol[var_num-1-i] * eq.pop
      end

      # quadratic terms
      var_num.times do |i|
        (var_num-i).times do |j|
          res += sol[var_num-1-i] * sol[var_num-1-i-j] * eq.pop
        end
      end

      abort "eq#{counter} is evaluated to 1!" if res.odd?
      counter += 1
    end
  end
end

puts 'solution correct!'
