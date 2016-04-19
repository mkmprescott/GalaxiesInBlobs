pro quicklook

cat = rsex('setest.cat')

select = where(cat.imaflags_iso ne 1)

cat = cat(select)

openw, 1, 'quicklook.reg'
for i=0, n_elements(cat)-1 do begin
printf, 1,'J2000; circle ', cat(i).alpha_j2000, cat(i).delta_j2000, ' 10p #text={', cat(i).number,'}'
endfor
close, 1
stop
end
