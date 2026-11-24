import fs from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";
import MarkdownIt from "markdown-it";
import { katex } from "@mdit/plugin-katex";
import GithubSlugger from "github-slugger";
import { parse as parseYaml } from "yaml";

export const root = fileURLToPath(new URL("../", import.meta.url));
export const repository = "https://github.com/lixiang90/vibemath-NumberTheory";
export const pages = "https://lixiang90.github.io/vibemath-NumberTheory";
export const escape = MarkdownIt().utils.escapeHtml;
const posix = path.posix;
const math = { delimiters: "all", mathFence: true, throwOnError: true };
const markdown = new MarkdownIt({ html: false }).use(katex, math);
const labels = new MarkdownIt("zero").enable(["backticks", "emphasis", "escape"]).use(katex, math);
export const inline = (value) => labels.renderInline(value);

export function encodePath(value) {
  return value.split("/").map(encodeURIComponent).join("/")
    .replace(/[!'()*]/g, c => "%" + c.charCodeAt(0).toString(16).toUpperCase());
}

export function href(from, to) {
  return encodePath(posix.relative(posix.dirname(from), to));
}

export function sourceUrl(source, directory = false) {
  return `${repository}/${directory ? "tree" : "blob"}/main/${encodePath(source)}`;
}

function plainText(token) {
  if (!token.children) return token.content;
  return token.children.map(plainText).join("");
}

function loadPage(source, output) {
  try {
    let body = fs.readFileSync(path.join(root, source), "utf8").replace(/^\uFEFF/, "");
    const front = /^---[ \t]*\r?\n([\s\S]*?)\r?\n---[ \t]*(?:\r?\n|$)/.exec(body);
    const meta = front ? parseYaml(front[1]) ?? {} : {};
    if (front) body = body.slice(front[0].length);
    const tokens = markdown.parse(body, {});
    const firstHeading = tokens.findIndex(token => token.type === "heading_open" && token.tag === "h1");
    const headingTitle = firstHeading >= 0 ? tokens[firstHeading + 1].content : "";
    const title = String(meta.title || headingTitle || posix.parse(source).name);
    return {
      source, output, title, body,
      author: [meta.author ?? []].flat().map(String).join("、"),
      date: String(meta.date ?? ""),
      abstract: String(meta.abstract ?? ""),
      example: meta.kind === "example",
      formula: String(meta.formula ?? "")
    };
  } catch (error) {
    throw new Error(`${source}: ${error.message}`, { cause: error });
  }
}

function* filesUnder(directory) {
  for (const entry of fs.readdirSync(path.join(root, directory), { withFileTypes: true })) {
    if (entry.name.startsWith(".")) continue;
    const source = posix.join(directory, entry.name);
    if (entry.isDirectory()) yield* filesUnder(source);
    else if (entry.isFile()) yield source;
  }
}

export function collect() {
  const topics = JSON.parse(fs.readFileSync(path.join(root, "site/topics.json"), "utf8"));
  const routes = new Map([
    [".", "index.html"], ["README.md", "index.html"],
    ["CONTRIBUTING.md", "guide/index.html"], ["templates/paper.md", "templates/paper.md"]
  ]);
  const assets = [];
  const guide = loadPage("CONTRIBUTING.md", "guide/index.html");
  const titles = new Intl.Collator("zh-Hans", { numeric: true });
  const outputs = new Set(routes.values());
  const reserved = new Set(["assets", "dist", "guide", "node_modules", "site", "templates"]);

  for (const topic of topics) {
    if (typeof topic.id !== "string" || !topic.id || topic.id.startsWith(".") ||
        /[\\/:*?"<>|]/.test(topic.id) || reserved.has(topic.id.toLowerCase())) {
      throw new Error(`Invalid or duplicate topic directory: ${topic.id}`);
    }
    reserved.add(topic.id.toLowerCase());
    topic.output = `${topic.id}/index.html`;
    topic.papers = [];
    routes.set(topic.id, topic.output);
    routes.set(`${topic.id}/README.md`, topic.output);
    if (!fs.existsSync(path.join(root, topic.id))) continue;
    for (const source of filesUnder(topic.id)) {
      if (source === `${topic.id}/README.md`) continue;
      const name = posix.relative(topic.id, source);
      const isPaper = /\.md$/i.test(name);
      const output = isPaper
        ? `${topic.id}/articles/${name.slice(0, -3)}.html`
        : `${topic.id}/files/${name}`;
      if (outputs.has(output)) throw new Error(`Duplicate output path: ${output}`);
      outputs.add(output);
      routes.set(source, output);
      if (isPaper) topic.papers.push(loadPage(source, output));
      else assets.push({ source, output });
    }
    topic.papers.sort((a, b) => b.date.localeCompare(a.date) || titles.compare(a.title, b.title) || a.source.localeCompare(b.source));
  }
  return { topics, routes, assets, guide };
}

// Resolve against the Markdown file before making the link relative to the generated page.
function rewrite(url, { page, routes }, attribute) {
  if (!url || /^(?:[a-z][a-z\d+.-]*:|\/\/|[?#])/i.test(url)) return url;
  const [, pathname, suffix] = /^([^?#]*)(.*)$/s.exec(url);
  const decoded = decodeURIComponent(pathname);
  const base = decoded.startsWith("/") ? "." : posix.dirname(page.source);
  const target = posix.normalize(posix.join(base, decoded)).replace(/^\/|\/$/g, "") || ".";
  if (target === ".." || target.startsWith("../")) throw new Error(`Link outside repository: ${url}`);
  const output = routes.get(target) ?? routes.get(posix.join(target, "README.md"));
  if (output) return href(page.output, output) + suffix;
  const stat = fs.statSync(path.join(root, target), { throwIfNoEntry: false });
  if (!stat) throw new Error(`Missing local link target: ${url}`);
  if (attribute === "src" && stat.isFile()) return `${repository}/raw/refs/heads/main/${encodePath(target)}${suffix}`;
  return sourceUrl(target, stat.isDirectory()) + suffix;
}

for (const [rule, attr] of [["link_open", "href"], ["image", "src"]]) {
  const previous = markdown.renderer.rules[rule] ?? ((tokens, i, options, env, self) => self.renderToken(tokens, i, options));
  markdown.renderer.rules[rule] = (tokens, i, options, env, self) => {
    const token = tokens[i];
    token.attrSet(attr, rewrite(token.attrGet(attr), env, attr));
    return previous(tokens, i, options, env, self);
  };
}

markdown.core.ruler.after("inline", "article_headings", ({ tokens, env }) => {
  if (!env.page) return;
  const slugger = new GithubSlugger();
  slugger.slug("main"); // Reserved for the skip link.
  env.titleId = slugger.slug(env.page.title);
  if (tokens[0]?.type === "heading_open" && tokens[0].tag === "h1" && tokens[1].content === env.page.title) {
    tokens.splice(0, 3);
  }
  for (const [i, token] of tokens.entries()) {
    if (token.type !== "heading_open") continue;
    const label = plainText(tokens[i + 1]);
    const id = slugger.slug(label);
    token.attrSet("id", id);
    if (token.tag === "h2") env.toc.push({ label, id });
  }
});

export function renderContent(page, routes) {
  try {
    const env = { page, routes, toc: [] };
    const body = markdown.render(page.body, env);
    return { titleId: env.titleId, toc: env.toc, body };
  } catch (error) {
    throw new Error(`${page.source}: ${error.message}`, { cause: error });
  }
}
