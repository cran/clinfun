c     calculates the Delong, Delong, Clarke-Pearson variance for area under ROC
c     nv, nn, nd are number of markers, normal and diseased subjects, n = nn+nd
      subroutine rocarea(n, nv, nn, nd, ratings, area, jkarea)
      integer n, nv, nn, nd
      double precision ratings(n,nv), area(nv), jkarea(n, nv)

      integer i, j, k
      double precision rn, rd, rnd, psi

      rn = dfloat((nn-1)*nd)
      rd = dfloat(nn*(nd-1))
      rnd = dfloat(nn*nd)
      do 50 k = 1, nv
         do 20 i = 1, nn
            do 10 j = nn+1, n
               psi = 0.0
               if (ratings(i,k) .lt. ratings(j,k)) psi = 1.0
               if (ratings(i,k) .eq. ratings(j,k)) psi = 0.5
c     calculate the sum of all psi(Xi, Xj) = I(Xi < Xj)
               area(k) = area(k) + psi
c     calculate the contribution of observation i (normal)
               jkarea(i,k) = jkarea(i,k) + psi
c     calculate the contribution of observation j (diseased)
               jkarea(j,k) = jkarea(j,k) + psi
 10         continue
 20      continue
c     calculate the area if observation i (normal) is left out
         do 30 i = 1, nn
            jkarea(i,k) = (area(k) - jkarea(i,k))/rn
 30      continue
c     calculate the area if observation j (diseased) is left out
         do 40 j = nn+1, n
            jkarea(j,k) = (area(k) - jkarea(j,k))/rd
 40      continue
c     calculate the area
         area(k) = area(k)/rnd
 50   continue

      return
      end
