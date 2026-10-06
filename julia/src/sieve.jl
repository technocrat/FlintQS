# SPDX-License-Identifier: GPL-2.0-or-later

function sieve!(sv::Vector{UInt8}, st::SIQSState)
    fill!(sv, 0)
    M = st.params.M
    len = length(sv)
    fb = st.fb
    @inbounds for ip in st.params.firstprime+1:length(fb.primes)
        p = fb.primes[ip]
        sz = fb.sizes[ip]
        a = st.r1[ip]
        b = st.r2[ip]
        for j in mod(a + M, p)+1:p:len
            sv[j] += sz
        end
        if b != a
            for j in mod(b + M, p)+1:p:len
                sv[j] += sz
            end
        end
    end
    return sv
end
