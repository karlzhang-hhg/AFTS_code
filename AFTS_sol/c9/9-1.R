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

da = read.table("data/m-fac-ex-9008.txt", header = T)
dim(da)
xmtx = cbind(rep(1, 228), da[, 14])
rtn = as.matrix(da[, 1:13])
xit.hat = solve(t(xmtx) %*% xmtx) %*% (t(xmtx) %*% rtn)
beta.hat = t(xit.hat[2, ])
E.hat = rtn - xmtx %*% xit.hat
D.hat = diag(t(E.hat) %*% E.hat / (228 - 2))
r.square = 1 - (228 - 2) * D.hat / diag(t(rtn) %*% rtn)
t(rbind(beta.hat, sqrt(D.hat), r.square))

