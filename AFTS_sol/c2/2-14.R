#
# Copyright (c) 2015 by Yuchao Zhao, Xiaoye Meng.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

require(TSA)

da = read.table("data/sp5may.dat", header = T)
y = diff(da[, 1]) * 100
x = diff(da[, 2]) * 100
cat("correlation =", cor(y, x), "\n")
m1 = lm(y ~ x)
summary(m1)
Box.test(m1$residuals, lag = 10, type = 'Ljung')
eacf(m1$residuals)
m2 = arima(y, order = c(1, 0, 1), xreg = x, include.mean = F)
m2
Box.test(m2$residuals, lag = 10, type = 'Ljung')

