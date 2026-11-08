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
          r.textContent = 'Requisition entered. You are on the waitlist in order of receipt. You will be written to when allocation opens. No payment has been taken.';
          f.reset();
        })
        .catch(function () {
          r.textContent = 'The requisition could not be filed. Write to ' + CONFIG.mailto + ' and it will be entered by hand.';
        })
        .then(function () { btn.disabled = false; });
    });
  });
})();
