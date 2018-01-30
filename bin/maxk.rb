#!/usr/bin/env ruby
# usage: find the max number of variables one can keep
# arguments:
#   1) var_num: number of variables
#   2) eq_num: number of eqs
#   3) deg: degree of Macaulay

abort 'invalid arguments' unless ARGV.size >= 3
var_num, eq_num, deg, _ = ARGV.map(&:to_i)

def binom2(n)
  n * (n-1) / 2
end

def binom3(n)
  n * (n-1) * (n-2) / 6
end

def binom4(n)
  n * (n-1) * (n-2) * (n-3) / 24
end

def nlnum(n, deg, k) 
  c = binom2(k) + binom3(k) + binom2(k) * (n-k)

  if 4 == deg
    c += binom4(k) + binom3(k) * (n-k) + binom2(k) * binom2(n-k)
  end

  c
end

mac_eq_num = ( (var_num + 1) + ((deg == 4) ? binom2(var_num) : 0) ) * eq_num
mac_term_num = 1 + var_num + binom2(var_num) + binom3(var_num)
mac_term_num += binom4(var_num) if deg == 4
mac_memsize = mac_eq_num.to_f * mac_term_num / 8 / (2 ** 30)

mac_indep_eq_num = mac_eq_num - ((3 == deg) ? 0 : eq_num * (eq_num+1) / 2)

# compute max k
max_k, nl_term_num = var_num.downto(0) do |k|
  ntnum = nlnum(var_num, deg, k)
  break k, ntnum if mac_indep_eq_num - ntnum >= k
end

puts <<-EOF
  number of variables: #{var_num}
  number of equations: #{eq_num}
  degree of Macaulay matrix: #{deg}
  number of rows in Macaulay matrix: #{mac_eq_num}
  number of cols in Macaulay matrix: #{mac_term_num}
  number of terms to eliminate: #{nl_term_num}
  max number of variables one can keep: #{max_k}
EOF
