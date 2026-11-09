/* Juris waitlist: requisition form JPI-R-1, source attribution, analytics.
   Configure once the Google Form and Cloudflare site exist; until then the form falls back to email. */
(function () {
  var CONFIG = {
    // Google Form "formResponse" URL, e.g. https://docs.google.com/forms/d/e/XXXX/formResponse
    formAction: '',
    // entry.NNNN ids from the form's pre-filled link, one per field below
    entries: { email: '', instrument: '', finish: '', fee: '', profession: '', heard: '', source: '' },
    // Cloudflare Web Analytics beacon token
    cfToken: '',
    mailto: 'info@jurisinstruments.com'
  };

  // ---- analytics
  if (CONFIG.cfToken) {
    var s = document.createElement('script'); s.defer = true;
    s.src = 'https://static.cloudflareinsights.com/beacon.min.js';
    s.setAttribute('data-cf-beacon', JSON.stringify({ token: CONFIG.cfToken }));
    document.head.appendChild(s);
  }

  // ---- source attribution: ?ref= / utm_source, else the referring site; kept for the visit
  function store(k, v) { try { sessionStorage.setItem(k, v); } catch (e) {} }
  function read(k) { try { return sessionStorage.getItem(k) || ''; } catch (e) { return ''; } }
  try {
    var q = new URLSearchParams(location.search);
    var ref = q.get('ref') || q.get('utm_source');
    if (ref) store('jpi_ref', ref + (q.get('utm_campaign') ? '/' + q.get('utm_campaign') : ''));
    else if (!read('jpi_ref') && document.referrer) {
      var h = new URL(document.referrer).hostname;
      if (h && h !== location.hostname) store('jpi_ref', h);
    }
  } catch (e) {}

  function entered() { try { return localStorage.getItem('jpi_req') === '1'; } catch (e) { return false; } }
  function markEntered() { try { localStorage.setItem('jpi_req', '1'); } catch (e) {} }

  // ---- forms
  document.querySelectorAll('form[data-waitlist]').forEach(function (f) {
    var r = f.querySelector('.receipt'), e = f.querySelector('input[type=email]'), btn = f.querySelector('button[type=submit]');
    if (entered()) { r.hidden = false; r.textContent = 'A requisition from this device is already on file. Another may be filed.'; }
    f.addEventListener('submit', function (ev) {
      ev.preventDefault();
      r.hidden = false;
      if (!e.value || !e.checkValidity()) { r.textContent = 'Field 1 is incomplete. A correspondence address is required.'; e.focus(); return; }
      var data = {};
      Array.prototype.forEach.call(f.elements, function (el) { if (el.name) data[el.name] = el.value; });
      data.source = read('jpi_ref') || 'direct';
      if (!CONFIG.formAction) {
        var body = Object.keys(data).map(function (k) { return k + ': ' + data[k]; }).join('\n');
        location.href = 'mailto:' + CONFIG.mailto + '?subject=' + encodeURIComponent('Requisition JPI-R-1') + '&body=' + encodeURIComponent(body);
        r.textContent = 'Your mail program has been asked to file this requisition. Send the message to complete it.';
        return;
      }
      var p = new URLSearchParams();
      Object.keys(CONFIG.entries).forEach(function (k) { if (CONFIG.entries[k] && data[k] != null) p.append(CONFIG.entries[k], data[k]); });
      btn.disabled = true; r.textContent = 'Filing…';
      fetch(CONFIG.formAction, { method: 'POST', mode: 'no-cors', body: p })
        .then(function () {
          markEntered();
          r.innerHTML = '<div style="width:220px;max-width:70%;margin:6px 0 14px;text-transform:none;letter-spacing:normal"><svg class="motto-lk" viewBox="0 0 148 64" role="img" aria-label="Juris Professional Instruments. Res ipsa loquitur."><g transform="translate(0,0)"><clipPath id="rc-c"><path d="M10,0 H22 L32,10 V22 L22,32 H10 L0,22 V10 Z"/></clipPath><g fill="currentColor" clip-path="url(#rc-c)"><rect x="0" y="0" width="9" height="9"/><rect x="0" y="11" width="9" height="10"/><rect x="0" y="23" width="9" height="9"/><rect x="11" y="0" width="10" height="9"/><rect x="11" y="11" width="10" height="10"/><rect x="11" y="23" width="10" height="9"/><rect x="23" y="0" width="9" height="9"/><rect x="23" y="11" width="9" height="10"/><rect x="23" y="23" width="9" height="9"/></g></g><mask id="rc-m" maskUnits="userSpaceOnUse" x="40" y="-2" width="112" height="36"><rect x="40" y="-2" width="112" height="36" fill="#fff"/><rect x="40" y="9" width="112" height="2" fill="#000"/><rect x="40" y="21" width="112" height="2" fill="#000"/></mask><text x="44" y="32" font-family="Barlow Condensed, Arial Narrow, sans-serif" font-weight="700" font-size="45.71" textLength="104" lengthAdjust="spacing" fill="currentColor" mask="url(#rc-m)">JURIS</text><text x="44.4" y="44" font-family="Barlow Condensed, Arial Narrow, sans-serif" font-weight="500" font-size="8.57" letter-spacing=".6" fill="currentColor">PROFESSIONAL</text><text x="148" y="44" text-anchor="end" font-family="Barlow Condensed, Arial Narrow, sans-serif" font-weight="500" font-size="8.57" letter-spacing=".6" fill="currentColor">INSTRUMENTS</text><rect x="0" y="50" width="148" height="0.8" fill="currentColor"/><text x="44" y="61.5" text-anchor="start" font-family="Old Standard TT, Century Schoolbook, Georgia, serif" font-style="italic" font-size="10.5" fill="currentColor">Res ipsa loquitur</text></svg></div>Requisition entered. You are on the waitlist in order of receipt. You will be written to when allocation opens. No payment has been taken.';
          f.reset();
        })
        .catch(function () {
          r.textContent = 'The requisition could not be filed. Write to ' + CONFIG.mailto + ' and it will be entered by hand.';
        })
        .then(function () { btn.disabled = false; });
    });
  });
})();
